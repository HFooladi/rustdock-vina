//! Monte Carlo search for molecular docking: AutoDock Vina's algorithm.
//!
//! The search itself, a port of Vina 1.2.7's, lives in the private
//! `vina_search` module; this module sets it up and turns its poses into
//! [`DockingResult`]s.

use nalgebra::Vector3;
use rand::prelude::*;

use crate::forcefield::ForceField;
use crate::grid::{pairwise_energy_and_gradient, GridSet};
use crate::molecule::Molecule;
use crate::optimization::vina_search::{
    vina_global_steps, vina_local_steps, ExactField, GridField, InterField, SearchBox,
    SearchSettings, Searcher,
};
use crate::optimization::{DockingResult, OptimizationError, Optimizer};

/// Parameters for the Monte Carlo search. Defaults are Vina's.
#[derive(Debug, Clone)]
pub struct MonteCarloParams {
    /// Independent Monte Carlo chains, run in parallel (Vina's
    /// `exhaustiveness`).
    pub exhaustiveness: usize,

    /// Metropolis temperature, in kcal/mol (Vina: 1.2).
    pub temperature: f64,

    /// Size of one mutation: Angstroms of translation, and the rotation that
    /// moves a typical atom about as far (Vina: 2).
    pub mutation_amplitude: f64,

    /// Cap that softens clashes while the search explores (Vina's
    /// `hunt_cap`: 10). Positive energies `e` count as `e * cap / (cap + e)`.
    pub hunt_cap: f64,

    /// Minimum heavy-atom RMSD between reported modes, in Angstroms
    /// (Vina's `--min_rmsd`: 1).
    pub min_rmsd: f64,

    /// BFGS steps per local optimization; `None` uses Vina's
    /// `(25 + atoms) / 3`.
    pub local_steps: Option<usize>,

    /// Whether each chain starts from a random pose in the box. Disable only
    /// to refine a pose you already have.
    pub randomize_initial_pose: bool,

    /// Optional master RNG seed for reproducible runs.
    ///
    /// `None` (the default) draws fresh entropy. `Some(seed)` makes the
    /// search fully deterministic: each chain's seed is mixed from the master
    /// seed and the chain index.
    pub seed: Option<u64>,

    /// Score ligand-receptor energy from precomputed grids during the search.
    ///
    /// Each pose then costs one interpolation per ligand atom instead of a sum
    /// over every nearby receptor atom. Reported energies are always rescored
    /// exactly. Ignored for force fields that cannot be gridded
    /// (see [`ForceField::grid_cutoff`]).
    pub use_grid: bool,

    /// Grid point spacing in Angstroms (Vina's default is 0.375).
    pub grid_spacing: f64,
}

impl Default for MonteCarloParams {
    fn default() -> Self {
        Self {
            exhaustiveness: 8,
            temperature: 1.2,
            mutation_amplitude: 2.0,
            hunt_cap: 10.0,
            min_rmsd: 1.0,
            local_steps: None,
            randomize_initial_pose: true,
            seed: None,
            use_grid: true,
            grid_spacing: 0.375,
        }
    }
}

/// Implementation of Monte Carlo search for molecular docking
#[derive(Debug, Clone, Default)]
pub struct MonteCarlo {
    pub params: MonteCarloParams,
}

impl MonteCarlo {
    /// Create a new Monte Carlo optimizer with default parameters
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a new Monte Carlo optimizer with custom parameters
    pub fn with_params(params: MonteCarloParams) -> Self {
        Self { params }
    }

    /// Vina's Monte Carlo step count per chain for `ligand`:
    /// `70 * 3 * (50 + atoms + 10 * (6 + torsions)) / 2`.
    pub fn default_steps(ligand: &Molecule) -> usize {
        vina_global_steps(ligand)
    }

    /// Dock `ligand`: up to `num_modes` distinct poses, best first, from
    /// `steps` Monte Carlo steps per chain (`0` for Vina's default).
    #[allow(clippy::too_many_arguments)]
    fn dock(
        &self,
        ligand: &Molecule,
        receptor: &Molecule,
        forcefield: &dyn ForceField,
        center: Vector3<f64>,
        box_size: Vector3<f64>,
        num_modes: usize,
        steps: usize,
    ) -> Result<Vec<DockingResult>, OptimizationError> {
        // Receptor atoms too far from the box can never interact with any
        // pose inside it, so drop them once.
        let pruned = prune_receptor(receptor, &center, &box_size, ligand);
        let search_box = SearchBox::new(&center, &box_size, self.params.grid_spacing);
        let grids = if self.params.use_grid {
            GridSet::build(
                ligand,
                &pruned,
                forcefield,
                &center,
                &box_size,
                0.0,
                self.params.grid_spacing,
            )?
        } else {
            None
        };
        let exact = ExactField::new(&pruned, forcefield, search_box);
        let grid_field = grids.as_ref().map(GridField);
        let field: &dyn InterField = match &grid_field {
            Some(g) => g,
            None => &exact,
        };

        let searcher = Searcher::new(ligand, forcefield, search_box)?;
        let settings = SearchSettings {
            chains: self.params.exhaustiveness.max(1),
            global_steps: if steps == 0 {
                vina_global_steps(ligand)
            } else {
                steps
            },
            local_steps: self
                .params
                .local_steps
                .unwrap_or_else(|| vina_local_steps(ligand)),
            num_saved: num_modes.max(1),
            min_rmsd: self.params.min_rmsd,
            temperature: self.params.temperature,
            mutation_amplitude: self.params.mutation_amplitude,
            hunt_cap: self.params.hunt_cap,
            randomize: self.params.randomize_initial_pose,
        };
        // Each chain's seed is mixed from (master, chain index), so the set
        // of chains is reproducible and no two master seeds share chains.
        let master = self.params.seed.unwrap_or_else(|| rand::thread_rng().gen());
        let poses = searcher.search(field, &exact, &settings, |i| run_seed(master, i as u64))?;

        // As Vina reports a mode: its intermolecular energy plus its internal
        // energy relative to the best mode's, which stands in for the
        // unbound ligand. The best mode's value is therefore its
        // intermolecular energy alone.
        let mut results = Vec::with_capacity(poses.len());
        let mut best_intra = None;
        for pose in &poses {
            let coords = searcher.coords_of(&pose.conf);
            let molecule = searcher.model.to_molecule(ligand, &pose.conf, &coords);
            let inter = intermolecular_energy(&molecule, &pruned, forcefield)?;
            let intra = searcher.intra_energy(&coords)?;
            let reference = *best_intra.get_or_insert(intra);
            results.push(DockingResult {
                energy: inter + intra - reference,
                molecule,
                energy_components: vec![
                    ("Intermolecular".to_string(), inter),
                    ("Intramolecular".to_string(), intra),
                ],
                rmsd: None,
            });
        }

        if let Some((first, rest)) = results.split_first_mut() {
            for result in rest {
                result.rmsd = Some(calculate_rmsd(&result.molecule, &first.molecule));
            }
        }
        Ok(results)
    }
}

/// How the search evaluates ligand-receptor energy: from precomputed grids
/// when available, otherwise by exact pairwise summation.
#[derive(Clone, Copy)]
pub(crate) struct InterModel<'a> {
    receptor: &'a Molecule,
    grids: Option<&'a GridSet>,
}

impl<'a> InterModel<'a> {
    /// Exact pairwise scoring against `receptor`.
    pub(crate) fn exact(receptor: &'a Molecule) -> Self {
        Self {
            receptor,
            grids: None,
        }
    }

    /// Grid scoring, with `receptor` behind it.
    #[cfg(test)]
    pub(crate) fn with_grids(receptor: &'a Molecule, grids: &'a GridSet) -> Self {
        Self {
            receptor,
            grids: Some(grids),
        }
    }

    pub(crate) fn energy(
        &self,
        ligand: &Molecule,
        forcefield: &dyn ForceField,
    ) -> Result<f64, OptimizationError> {
        match self.grids {
            Some(grids) => Ok(grids.energy(ligand, forcefield)?),
            None => intermolecular_energy(ligand, self.receptor, forcefield),
        }
    }

    /// As [`energy`](Self::energy), additionally adding each ligand atom's
    /// gradient `dE/dr_i` into `grads`.
    pub(crate) fn energy_and_gradients(
        &self,
        ligand: &Molecule,
        forcefield: &dyn ForceField,
        grads: &mut [Vector3<f64>],
    ) -> Result<f64, OptimizationError> {
        if let Some(grids) = self.grids {
            return Ok(grids.energy_and_gradients(ligand, forcefield, grads)?);
        }
        let mut total = 0.0;
        for (atom, grad) in ligand.atoms.iter().zip(grads) {
            let (energy, gradient) =
                pairwise_energy_and_gradient(atom, &self.receptor.atoms, forcefield)?;
            total += energy;
            *grad += gradient;
        }
        Ok(total)
    }
}

/// Intermolecular (ligand-receptor) interaction energy.
///
/// This is the quantity the reported affinity is derived from.
pub(crate) fn intermolecular_energy(
    ligand: &Molecule,
    receptor: &Molecule,
    forcefield: &dyn ForceField,
) -> Result<f64, OptimizationError> {
    let mut total = 0.0;
    for lig_atom in &ligand.atoms {
        for rec_atom in &receptor.atoms {
            let distance = lig_atom.distance(rec_atom);
            if distance > 0.01 {
                total += forcefield.atom_pair_energy(lig_atom, rec_atom, distance)?;
            }
        }
    }
    Ok(total)
}

/// Intramolecular (ligand self-interaction) energy over the given pairs.
///
/// Without this, rotating a torsion is free: the search can fold the ligand
/// through itself and pay nothing, because only receptor contacts are scored.
/// It matters exactly when torsions are being sampled, which is always.
pub(crate) fn intramolecular_energy(
    ligand: &Molecule,
    pairs: &[(usize, usize)],
    forcefield: &dyn ForceField,
) -> Result<f64, OptimizationError> {
    let mut total = 0.0;
    for &(i, j) in pairs {
        let (a, b) = (&ligand.atoms[i], &ligand.atoms[j]);
        let distance = a.distance(b);
        if distance > 0.01 {
            total += forcefield.atom_pair_energy(a, b, distance)?;
        }
    }
    Ok(total)
}

/// As [`intramolecular_energy`], additionally adding each atom's gradient
/// into `grads`.
pub(crate) fn intramolecular_energy_and_gradients(
    ligand: &Molecule,
    pairs: &[(usize, usize)],
    forcefield: &dyn ForceField,
    grads: &mut [Vector3<f64>],
) -> Result<f64, OptimizationError> {
    let mut total = 0.0;
    for &(i, j) in pairs {
        let (a, b) = (&ligand.atoms[i], &ligand.atoms[j]);
        let delta = a.coordinates - b.coordinates;
        let distance = delta.norm();
        if distance > 0.01 {
            total += forcefield.atom_pair_energy(a, b, distance)?;
            let g = delta * (forcefield.atom_pair_derivative(a, b, distance)? / distance);
            grads[i] += g;
            grads[j] -= g;
        }
    }
    Ok(total)
}

impl Optimizer for MonteCarlo {
    /// The single best pose from a full search; `max_steps` is Monte Carlo
    /// steps per chain (`0` for Vina's default).
    fn optimize(
        &self,
        ligand: &Molecule,
        receptor: &Molecule,
        forcefield: &dyn ForceField,
        center: Vector3<f64>,
        box_size: Vector3<f64>,
        max_steps: usize,
    ) -> Result<DockingResult, OptimizationError> {
        self.dock(ligand, receptor, forcefield, center, box_size, 1, max_steps)?
            .into_iter()
            .next()
            .ok_or(OptimizationError::FailedToConverge)
    }

    /// Up to `num_poses` distinct binding modes, best first; `max_steps` is
    /// Monte Carlo steps per chain (`0` for Vina's default).
    fn generate_poses(
        &self,
        ligand: &Molecule,
        receptor: &Molecule,
        forcefield: &dyn ForceField,
        center: Vector3<f64>,
        box_size: Vector3<f64>,
        num_poses: usize,
        max_steps: usize,
    ) -> Result<Vec<DockingResult>, OptimizationError> {
        self.dock(
            ligand, receptor, forcefield, center, box_size, num_poses, max_steps,
        )
    }
}

/// Drop receptor atoms that cannot reach any ligand pose placed in the box.
///
/// A pose's centroid is confined to the box, so no ligand atom can be further
/// than (half-box + ligand radius) from the box centre; adding the interaction
/// cutoff gives a conservative bound. Grids built from the
/// pruned receptor are unaffected, since pruned atoms never reach the grid.
fn prune_receptor(
    receptor: &Molecule,
    center: &Vector3<f64>,
    box_size: &Vector3<f64>,
    ligand: &Molecule,
) -> Molecule {
    let half = box_size * 0.5;
    let margin = ligand_radius(ligand) + crate::forcefield::vina::CUTOFF;

    let mut pruned = receptor.clone();
    pruned.atoms.retain(|atom| {
        let d = atom.coordinates - center;
        d.x.abs() <= half.x + margin && d.y.abs() <= half.y + margin && d.z.abs() <= half.z + margin
    });

    // Bonds and torsions index into the original atom list, so they would be
    // wrong after pruning. Nothing in the scoring path uses receptor
    // connectivity — the xs typing flags were already baked into each atom at
    // parse time — so they are simply dropped.
    pruned.bonds.clear();
    pruned.torsions.clear();
    pruned
}

/// Seed for run `run` of a docking job with master seed `master`.
///
/// Simply adding the index made master `m` run `i + 1` identical to master
/// `m + 1` run `i`, so docks with nearby seeds silently repeated most of their
/// runs. Mixing both inputs through splitmix64 keeps every (master, run) pair
/// distinct.
fn run_seed(master: u64, run: u64) -> u64 {
    const GAMMA: u64 = 0x9E37_79B9_7F4A_7C15;
    splitmix64(splitmix64(master).wrapping_add(run.wrapping_add(1).wrapping_mul(GAMMA)))
}

/// The splitmix64 finalizer: a fast, well-mixed 64-bit bijection.
fn splitmix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Distance from the ligand's centroid to its furthest atom.
fn ligand_radius(ligand: &Molecule) -> f64 {
    ligand
        .center()
        .map(|c| {
            ligand
                .atoms
                .iter()
                .map(|a| (a.coordinates - c).norm())
                .fold(0.0, f64::max)
        })
        .unwrap_or(0.0)
}

/// Calculate root mean square deviation between two molecules
fn calculate_rmsd(mol1: &Molecule, mol2: &Molecule) -> f64 {
    if mol1.atoms.len() != mol2.atoms.len() {
        return f64::MAX; // Different number of atoms
    }

    let mut sum_sq_diff = 0.0;

    for i in 0..mol1.atoms.len() {
        let pos1 = &mol1.atoms[i].coordinates;
        let pos2 = &mol2.atoms[i].coordinates;

        let diff = pos1 - pos2;
        sum_sq_diff += diff.norm_squared();
    }

    (sum_sq_diff / mol1.atoms.len() as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    /// Nearby master seeds must not share runs. With `master + run`, masters
    /// 0..16 over 8 runs produced only 23 distinct seeds out of 128.
    #[test]
    fn test_run_seeds_do_not_collide_across_masters() {
        let seeds: HashSet<u64> = (0..16)
            .flat_map(|master| (0..8).map(move |run| run_seed(master, run)))
            .collect();
        assert_eq!(seeds.len(), 16 * 8);
    }
}
