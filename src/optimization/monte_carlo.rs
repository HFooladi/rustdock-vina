//! Monte Carlo search for molecular docking: AutoDock Vina's algorithm.
//!
//! The search itself, a port of Vina 1.2.7's, lives in the private
//! `vina_search` module; this module sets it up and turns its poses into
//! [`DockingResult`]s.

use nalgebra::Vector3;
use rand::prelude::*;

use crate::atom::{Atom, Element};
use crate::forcefield::ForceField;
use crate::grid::{pairwise_energy_and_gradient, GridSet, ReceptorGrids};
use crate::molecule::{FlexResidues, Molecule};
use crate::optimization::vina_search::{
    vina_global_steps, vina_local_steps, Caps, ExactField, GridField, InterField, SearchBox,
    SearchSettings, Searcher,
};
use crate::optimization::{DockingResult, OptimizationError, Optimizer};
use crate::scoring::score_pose_with_flex;

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

    /// Cap that softens the ligand's internal clashes while the search
    /// explores (Vina's `hunt_cap[0]` and `[2]`: 10). Positive energies `e`
    /// count as `e * cap / (cap + e)`.
    pub hunt_cap: f64,

    /// The same cap on each ligand atom's energy against the receptor
    /// (Vina's `hunt_cap[1]`: 10). Vina's `monte_carlo` defaults it to 1.5,
    /// but `Vina::global_search` sets all three caps to 10 before docking.
    pub hunt_cap_receptor: f64,

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
            hunt_cap_receptor: 10.0,
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
        vina_global_steps(ligand.atoms.len(), ligand.torsions.len())
    }

    /// The receptor and search box for docking with these parameters, with
    /// grids shared by every ligand docked into it (see [`Target`]).
    pub fn target(
        &self,
        receptor: Molecule,
        forcefield: &dyn ForceField,
        center: Vector3<f64>,
        box_size: Vector3<f64>,
    ) -> Result<Target, OptimizationError> {
        let grids = if self.params.use_grid {
            ReceptorGrids::new(
                &receptor,
                forcefield,
                &center,
                &box_size,
                0.0,
                self.params.grid_spacing,
            )?
        } else {
            None
        };
        Ok(Target {
            receptor,
            flex: None,
            center,
            box_size,
            grids,
        })
    }

    /// As [`target`](Self::target), with flexible side chains that every
    /// ligand docked into it moves along with itself. `receptor` is the
    /// rigid part, which must leave the side chains out.
    pub fn target_with_flex(
        &self,
        receptor: Molecule,
        flex: FlexResidues,
        forcefield: &dyn ForceField,
        center: Vector3<f64>,
        box_size: Vector3<f64>,
    ) -> Result<Target, OptimizationError> {
        let mut target = self.target(receptor, forcefield, center, box_size)?;
        target.flex = Some(flex);
        Ok(target)
    }

    /// Dock `ligand` into `target`, which must have been made by
    /// [`target`](Self::target) with the same force field: up to `num_modes`
    /// distinct poses, best first, from `steps` Monte Carlo steps per chain
    /// (`0` for Vina's default). The poses are the same as
    /// [`Optimizer::generate_poses`] gives.
    pub fn dock_into(
        &self,
        target: &Target,
        ligand: &Molecule,
        forcefield: &dyn ForceField,
        num_modes: usize,
        steps: usize,
    ) -> Result<Vec<DockingResult>, OptimizationError> {
        let (center, box_size) = (target.center, target.box_size);
        // Receptor atoms too far from the box can never interact with any
        // pose inside it, so drop them once.
        let pruned = prune_receptor(&target.receptor, &center, &box_size, ligand);
        let search_box = SearchBox::new(&center, &box_size, self.params.grid_spacing);
        let flex = target.flex.as_ref();
        let searcher = Searcher::new(ligand, flex, forcefield, search_box)?;
        let grids = match &target.grids {
            Some(grids) => {
                let mut movable = Molecule::new("movable");
                movable.atoms = searcher.movable_atoms().to_vec();
                Some(grids.grid_set(&movable, forcefield)?)
            }
            None => None,
        };
        let exact = ExactField::new(&pruned, forcefield, search_box);
        let grid_field = grids.as_ref().map(GridField);
        let field: &dyn InterField = match &grid_field {
            Some(g) => g,
            None => &exact,
        };

        let movable = searcher.model.n_movable();
        let torsions = searcher.model.n_torsions() + searcher.model.n_flex_torsions();
        let settings = SearchSettings {
            chains: self.params.exhaustiveness.max(1),
            global_steps: if steps == 0 {
                vina_global_steps(movable, torsions)
            } else {
                steps
            },
            local_steps: self
                .params
                .local_steps
                .unwrap_or_else(|| vina_local_steps(movable)),
            num_saved: num_modes.max(1),
            min_rmsd: self.params.min_rmsd,
            temperature: self.params.temperature,
            mutation_amplitude: self.params.mutation_amplitude,
            hunt_caps: Caps {
                receptor: self.params.hunt_cap_receptor,
                intra: self.params.hunt_cap,
            },
            randomize: self.params.randomize_initial_pose,
        };
        // Each chain's seed is mixed from (master, chain index), so the set
        // of chains is reproducible and no two master seeds share chains.
        let master = self.params.seed.unwrap_or_else(|| rand::thread_rng().gen());
        let poses = searcher.search(field, &exact, &settings, |i| run_seed(master, i as u64))?;

        // As Vina reports a mode: its intermolecular energy plus its internal
        // energy relative to the best mode's, which stands in for the
        // unbound ligand. The best mode's value is therefore its
        // intermolecular energy alone. With flexible side chains, their
        // energy with the ligand is intermolecular, and their energy with
        // the rigid receptor and with each other is internal.
        let mut results = Vec::with_capacity(poses.len());
        let mut best_intra = None;
        // Vina measures RMSD over the movable atoms: the ligand, then any
        // flexible side chain, in the search's own order.
        let movable_atoms = searcher.movable_atoms().to_vec();
        let n_movable_atoms = movable_atoms.len();
        let mut movable_coords: Vec<Vec<Vector3<f64>>> = Vec::with_capacity(poses.len());
        for pose in &poses {
            let coords = searcher.coords_of(&pose.conf);
            movable_coords.push(coords[..n_movable_atoms].to_vec());
            let molecule = searcher.model.to_molecule(ligand, &pose.conf, &coords);
            let mut inter = intermolecular_energy(&molecule, &pruned, forcefield)?;
            let mut intra = searcher.intra_energy(&coords)?;
            let moved = flex.map(|f| searcher.model.to_flex(f, &pose.conf, &coords));
            if let Some(moved) = &moved {
                let score = score_pose_with_flex(&molecule, &pruned, Some(moved), forcefield)
                    .map_err(|e| OptimizationError::Unknown(e.to_string()))?;
                inter += score.ligand_flex;
                intra += score.flex_receptor + score.flex_flex;
            }
            let reference = *best_intra.get_or_insert(intra);
            results.push(DockingResult {
                energy: inter + intra - reference,
                molecule,
                flex: moved.map(|f| f.molecule),
                energy_components: vec![
                    ("Intermolecular".to_string(), inter),
                    ("Intramolecular".to_string(), intra),
                    ("Unbound".to_string(), reference),
                ],
                rmsd_lb: None,
                rmsd_ub: None,
            });
        }

        // Each mode's distance from the best one, as Vina reports it.
        if let Some((first, rest)) = movable_coords.split_first() {
            for (result, coords) in results.iter_mut().skip(1).zip(rest) {
                result.rmsd_lb = Some(rmsd_lower_bound(&movable_atoms, coords, first));
                result.rmsd_ub = Some(rmsd_upper_bound(&movable_atoms, coords, first));
            }
        }
        Ok(results)
    }
}

/// A receptor and search box ready for docking, made by
/// [`MonteCarlo::target`].
///
/// Its affinity grids are built as ligands need them and shared by every
/// ligand docked into it, from any number of threads at once.
#[derive(Debug)]
pub struct Target {
    receptor: Molecule,
    flex: Option<FlexResidues>,
    center: Vector3<f64>,
    box_size: Vector3<f64>,
    grids: Option<ReceptorGrids>,
}

impl Target {
    /// The receptor, or its rigid part when there are flexible side
    /// chains.
    pub fn receptor(&self) -> &Molecule {
        &self.receptor
    }

    /// The flexible side chains, in their input conformation.
    pub fn flex(&self) -> Option<&FlexResidues> {
        self.flex.as_ref()
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
        let target = self.target(receptor.clone(), forcefield, center, box_size)?;
        self.dock_into(&target, ligand, forcefield, 1, max_steps)?
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
        let target = self.target(receptor.clone(), forcefield, center, box_size)?;
        self.dock_into(&target, ligand, forcefield, num_poses, max_steps)
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
/// Vina's two RMSD measures between two poses of the same molecule, as the
/// `rmsd l.b.` and `rmsd u.b.` columns of a docking run report them.
///
/// The upper bound pairs atom `i` with atom `i`; the lower bound pairs each
/// heavy atom with the nearest one of the same element, in whichever
/// direction gives the larger value, so a symmetric group that has merely
/// rotated onto itself does not count as movement. Hydrogens are left out
/// of both. Returns `None` if the two poses do not hold the same atoms.
///
/// With flexible side chains, Vina measures over the side chains too; this
/// compares the molecules it is given.
pub fn rmsd_bounds(a: &Molecule, b: &Molecule) -> Option<(f64, f64)> {
    if a.atoms.len() != b.atoms.len() {
        return None;
    }
    let coords = |m: &Molecule| -> Vec<Vector3<f64>> {
        m.atoms.iter().map(|atom| atom.coordinates).collect()
    };
    let (x, y) = (coords(a), coords(b));
    Some((
        rmsd_lower_bound(&a.atoms, &x, &y),
        rmsd_upper_bound(&a.atoms, &x, &y),
    ))
}

/// Vina's `rmsd_upper_bound`: heavy atoms paired by index.
///
/// `atoms` types both pose's coordinates, which are the search's movable
/// atoms in one order: the ligand, then any flexible side chain.
fn rmsd_upper_bound(atoms: &[Atom], a: &[Vector3<f64>], b: &[Vector3<f64>]) -> f64 {
    let mut sum = 0.0;
    let mut counted = 0usize;
    for (i, atom) in atoms.iter().enumerate() {
        if atom.atom_type.element() == Element::Hydrogen {
            continue;
        }
        sum += (a[i] - b[i]).norm_squared();
        counted += 1;
    }
    if counted == 0 {
        0.0
    } else {
        (sum / counted as f64).sqrt()
    }
}

/// One direction of Vina's `rmsd_lower_bound`: every heavy atom of `a`
/// against the nearest heavy atom of `b` of the same element.
fn rmsd_lower_bound_asymmetric(atoms: &[Atom], a: &[Vector3<f64>], b: &[Vector3<f64>]) -> f64 {
    let mut sum = 0.0;
    let mut counted = 0usize;
    for (i, atom) in atoms.iter().enumerate() {
        let element = atom.atom_type.element();
        if element == Element::Hydrogen {
            continue;
        }
        let nearest = atoms
            .iter()
            .enumerate()
            .filter(|(_, other)| {
                let other = other.atom_type.element();
                other == element && other != Element::Hydrogen
            })
            .map(|(j, _)| (a[i] - b[j]).norm_squared())
            .fold(f64::INFINITY, f64::min);
        // Unreachable while both poses hold the same atoms: `i` itself
        // always matches. Vina asserts the same.
        if nearest.is_finite() {
            sum += nearest;
            counted += 1;
        }
    }
    if counted == 0 {
        0.0
    } else {
        (sum / counted as f64).sqrt()
    }
}

/// Vina's `rmsd_lower_bound`: the larger of the two directions, so that a
/// pose cannot look close merely by crowding onto part of the other.
///
/// Matching each atom to the nearest of its element makes this blind to a
/// symmetry that swaps like atoms, where the upper bound counts the swap
/// as movement. Vina reports both.
fn rmsd_lower_bound(atoms: &[Atom], a: &[Vector3<f64>], b: &[Vector3<f64>]) -> f64 {
    rmsd_lower_bound_asymmetric(atoms, a, b).max(rmsd_lower_bound_asymmetric(atoms, b, a))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::AtomType;
    use std::collections::HashSet;

    fn atom_at(atom_type: AtomType, x: f64, y: f64, z: f64) -> Atom {
        Atom::new(
            atom_type,
            Vector3::new(x, y, z),
            "A".to_string(),
            1,
            "LIG".to_string(),
            1,
            'A',
            0.0,
        )
    }

    /// Two carbons that have swapped places are the same pose to the lower
    /// bound and a move to the upper one. This is the whole reason Vina
    /// reports both, so the two must not be computed alike.
    #[test]
    fn swapping_like_atoms_separates_the_two_bounds() {
        let atoms = [
            atom_at(AtomType::Carbon, 0.0, 0.0, 0.0),
            atom_at(AtomType::Carbon, 1.0, 0.0, 0.0),
        ];
        let a = vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)];
        let b = vec![Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0)];

        assert_eq!(rmsd_lower_bound(&atoms, &a, &b), 0.0);
        assert_eq!(rmsd_upper_bound(&atoms, &a, &b), 1.0);
    }

    /// Different elements cannot stand in for each other, so an oxygen
    /// moving is a move by both measures.
    #[test]
    fn unlike_atoms_do_not_match_each_other() {
        let atoms = [
            atom_at(AtomType::Carbon, 0.0, 0.0, 0.0),
            atom_at(AtomType::Oxygen, 1.0, 0.0, 0.0),
        ];
        let a = vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)];
        let b = vec![Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0)];

        assert_eq!(rmsd_lower_bound(&atoms, &a, &b), 1.0);
        assert_eq!(rmsd_upper_bound(&atoms, &a, &b), 1.0);
    }

    /// Vina measures over heavy atoms only; a hydrogen may move freely.
    #[test]
    fn hydrogens_are_left_out() {
        let atoms = [
            atom_at(AtomType::Carbon, 0.0, 0.0, 0.0),
            atom_at(AtomType::HydrogenD, 1.0, 0.0, 0.0),
        ];
        let a = vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)];
        let b = vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(9.0, 0.0, 0.0)];

        assert_eq!(rmsd_upper_bound(&atoms, &a, &b), 0.0);
        assert_eq!(rmsd_lower_bound(&atoms, &a, &b), 0.0);
    }

    /// The lower bound takes the worse direction, so a pose cannot look
    /// close by crowding every atom onto one of the other's.
    #[test]
    fn the_lower_bound_takes_the_worse_direction() {
        let atoms = [
            atom_at(AtomType::Carbon, 0.0, 0.0, 0.0),
            atom_at(AtomType::Carbon, 0.0, 0.0, 0.0),
        ];
        // Both of `b`'s atoms sit on the first of `a`'s.
        let a = vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(4.0, 0.0, 0.0)];
        let b = vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0)];

        // a -> b: the second atom is 4 away; b -> a: both find an atom at 0.
        assert_approx_eq::assert_approx_eq!(
            rmsd_lower_bound_asymmetric(&atoms, &b, &a),
            0.0,
            1e-12
        );
        assert_approx_eq::assert_approx_eq!(
            rmsd_lower_bound(&atoms, &a, &b),
            (16.0f64 / 2.0).sqrt(),
            1e-12
        );
    }

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
