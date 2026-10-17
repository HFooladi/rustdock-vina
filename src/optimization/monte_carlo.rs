//! Monte Carlo search algorithm for molecular docking

use nalgebra::{Unit, UnitQuaternion, Vector3};
use rand::prelude::*;
use rand::rngs::StdRng;
use rayon::prelude::*;
use std::f64::consts::PI;

use crate::forcefield::ForceField;
use crate::grid::{pairwise_energy_and_gradient, GridSet};
use crate::molecule::Molecule;
use crate::optimization::local::{LocalOptimizer, LocalOptimizerParams};
use crate::optimization::{DockingResult, OptimizationError, Optimizer};

/// Parameters for Monte Carlo search
#[derive(Debug, Clone)]
pub struct MonteCarloParams {
    /// Temperature parameter for Metropolis criterion
    pub temperature: f64,

    /// Probability of performing a translation move
    pub prob_translation: f64,

    /// Probability of performing a rotation move
    pub prob_rotation: f64,

    /// Probability of performing a torsion move
    pub prob_torsion: f64,

    /// Maximum translation step (in Angstroms)
    pub max_translation: f64,

    /// Maximum rotation step (in radians)
    pub max_rotation: f64,

    /// Maximum torsion step (in radians)
    pub max_torsion: f64,

    /// Number of steps without improvement before stopping
    pub steps_without_improvement: usize,

    /// Whether to use local optimization (BFGS) to refine poses
    pub use_local_optimization: bool,

    /// Whether each run starts from a random placement inside the search box.
    ///
    /// This is what makes the search global. Disable it only to refine a pose
    /// you already have (the equivalent of Vina's `--local_only`); with it off,
    /// docking a ligand that starts far from the site will not find the site.
    pub randomize_initial_pose: bool,

    /// Iterations of local minimization applied to each Monte Carlo proposal
    /// before the Metropolis test.
    ///
    /// Kept short because it runs on every step; the winning pose gets a
    /// full-length minimization at the end.
    pub refine_iterations: usize,

    /// Optional master RNG seed for reproducible runs.
    ///
    /// `None` (the default) draws fresh entropy from the OS, so each run
    /// produces different poses. `Some(seed)` makes [`MonteCarlo::optimize`]
    /// fully deterministic; [`MonteCarlo::generate_poses`] derives a distinct
    /// per-pose seed from the master seed so the parallel pose set is also
    /// reproducible.
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
            temperature: 1.2, // Vina's value; 0.3 froze the walk into its first basin
            prob_translation: 0.5,
            prob_rotation: 0.3,
            prob_torsion: 0.2,
            max_translation: 2.0,
            max_rotation: PI / 3.0, // 60 degrees
            max_torsion: PI / 3.0,  // 60 degrees
            steps_without_improvement: 1000,
            use_local_optimization: true, // Enable local optimization by default
            randomize_initial_pose: true,
            refine_iterations: 12,
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

    /// Generate a random translation within the search box
    fn random_translation(&self, rng: &mut StdRng) -> Vector3<f64> {
        let dx = (rng.gen::<f64>() - 0.5) * 2.0 * self.params.max_translation;
        let dy = (rng.gen::<f64>() - 0.5) * 2.0 * self.params.max_translation;
        let dz = (rng.gen::<f64>() - 0.5) * 2.0 * self.params.max_translation;

        Vector3::new(dx, dy, dz)
    }

    /// Generate a random rotation of bounded magnitude about a uniformly
    /// distributed axis.
    fn random_rotation(&self, rng: &mut StdRng) -> UnitQuaternion<f64> {
        let angle = (rng.gen::<f64>() - 0.5) * 2.0 * self.params.max_rotation;
        UnitQuaternion::from_axis_angle(&random_axis(rng), angle)
    }

    /// Place the ligand at a uniformly random position and orientation in the
    /// box, with all torsions randomized.
    ///
    /// Placement targets the centroid, so this always succeeds — no rejection
    /// loop and no fallback to the input pose, which would quietly reintroduce
    /// a dependence on where the caller happened to put the ligand.
    fn random_pose_in_box(
        &self,
        ligand: &Molecule,
        center: &Vector3<f64>,
        box_size: &Vector3<f64>,
        rng: &mut StdRng,
    ) -> Molecule {
        let mut candidate = ligand.clone();

        for i in 0..candidate.torsions.len() {
            let angle = (rng.gen::<f64>() - 0.5) * 2.0 * PI;
            let _ = candidate.rotate_torsion(i, angle);
        }

        candidate.rotate_about_center(&random_orientation(rng));

        let Ok(current_center) = candidate.center() else {
            return ligand.clone();
        };
        let half = box_size * 0.5;
        let target = Vector3::new(
            center.x + (rng.gen::<f64>() - 0.5) * 2.0 * half.x,
            center.y + (rng.gen::<f64>() - 0.5) * 2.0 * half.y,
            center.z + (rng.gen::<f64>() - 0.5) * 2.0 * half.z,
        );
        candidate.translate(&(target - current_center));

        candidate
    }

    /// The quantity the search minimizes: receptor interaction plus the
    /// ligand's own internal energy, so that torsional moves cannot buy
    /// receptor contacts by folding the ligand into itself.
    fn objective(
        &self,
        ligand: &Molecule,
        inter: InterModel,
        forcefield: &dyn ForceField,
        intra_pairs: &[(usize, usize)],
    ) -> Result<f64, OptimizationError> {
        Ok(inter.energy(ligand, forcefield)?
            + intramolecular_energy(ligand, intra_pairs, forcefield)?)
    }

    /// Build the search grids for this ligand and box, if enabled and the
    /// force field allows it.
    ///
    /// The grid spans the box widened by the ligand's radius, the furthest an
    /// atom can sit from a centroid confined to the box (barring torsional
    /// extension, which the exact fallback in [`GridSet`] covers).
    fn build_grids(
        &self,
        ligand: &Molecule,
        receptor: &Molecule,
        forcefield: &dyn ForceField,
        center: &Vector3<f64>,
        box_size: &Vector3<f64>,
    ) -> Result<Option<GridSet>, OptimizationError> {
        if !self.params.use_grid {
            return Ok(None);
        }
        Ok(GridSet::build(
            ligand,
            receptor,
            forcefield,
            center,
            box_size,
            ligand_radius(ligand),
            self.params.grid_spacing,
        )?)
    }

    /// One Monte Carlo run, searching on `grids` when given and on exact
    /// pairwise energy otherwise.
    #[allow(clippy::too_many_arguments)]
    fn run(
        &self,
        ligand: &Molecule,
        receptor: &Molecule,
        grids: Option<&GridSet>,
        forcefield: &dyn ForceField,
        center: Vector3<f64>,
        box_size: Vector3<f64>,
        max_steps: usize,
    ) -> Result<DockingResult, OptimizationError> {
        let inter = InterModel { receptor, grids };
        let mut rng = match self.params.seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };
        // Start each run from a random placement inside the box rather than
        // from the input coordinates. Without this every replica explores the
        // same basin, "exhaustiveness" buys almost nothing, and a ligand handed
        // in at its answer is simply handed back.
        let mut current_molecule = if self.params.randomize_initial_pose {
            self.random_pose_in_box(ligand, &center, &box_size, &mut rng)
        } else {
            ligand.clone()
        };

        // Topology-derived, so computed once and reused for every evaluation.
        let intra_pairs = ligand.intramolecular_pairs();

        // Calculate initial energy
        let mut current_energy =
            self.objective(&current_molecule, inter, forcefield, &intra_pairs)?;
        let mut best_energy = current_energy;
        let mut best_molecule = current_molecule.clone();

        let mut steps_since_improvement = 0;

        // Vina's basin-hopping: perturb, minimize, then apply the Metropolis
        // test to the *minimized* energy. Testing the raw perturbed energy
        // instead makes the walk hop between unrelaxed configurations and
        // wastes most proposals, because almost any random nudge raises the
        // energy before relaxation has a chance to recover it.
        let refiner = LocalOptimizer::with_params(LocalOptimizerParams {
            max_iterations: self.params.refine_iterations,
            ..LocalOptimizerParams::default()
        });
        let bounds = Some((&center, &box_size));

        for _step in 0..max_steps {
            // Create a copy of the current state
            let mut new_molecule = current_molecule.clone();

            // Perturb every degree of freedom at once, as Vina does, rather
            // than picking a single one per step.
            new_molecule.translate(&self.random_translation(&mut rng));
            new_molecule.rotate_about_center(&self.random_rotation(&mut rng));
            for i in 0..new_molecule.torsions.len() {
                let delta = (rng.gen::<f64>() - 0.5) * 2.0 * self.params.max_torsion;
                let _ = new_molecule.rotate_torsion(i, delta);
            }

            // Check if the molecule is still within the search box
            if !new_molecule.center_within_box(&center, &box_size) {
                continue;
            }

            let new_energy = if self.params.use_local_optimization {
                refiner
                    .minimize_with_model(&mut new_molecule, inter, forcefield, bounds, &intra_pairs)
                    .unwrap_or(f64::INFINITY)
            } else {
                self.objective(&new_molecule, inter, forcefield, &intra_pairs)?
            };

            // Decide whether to accept the move (Metropolis criterion)
            let accept = if new_energy <= current_energy {
                true
            } else {
                let boltzmann_factor =
                    ((current_energy - new_energy) / self.params.temperature).exp();
                rng.gen::<f64>() < boltzmann_factor
            };

            if accept {
                current_molecule = new_molecule;
                current_energy = new_energy;

                // Check if this is a new best solution
                if new_energy < best_energy {
                    best_energy = new_energy;
                    best_molecule = current_molecule.clone();
                    steps_since_improvement = 0;
                } else {
                    steps_since_improvement += 1;
                }
            } else {
                steps_since_improvement += 1;
            }

            // Check termination criteria
            if steps_since_improvement >= self.params.steps_without_improvement {
                break;
            }
        }

        // Polish the winner with a full-length minimization; the in-loop one is
        // deliberately truncated for speed. Like Vina's final refinement, this
        // runs on exact energies: grid interpolation undersamples the narrow
        // attractive wells, which matters most in exactly the tight contacts a
        // good pose has.
        if self.params.use_local_optimization {
            let local_optimizer = LocalOptimizer::new();
            // Only the pose is kept; its energies are rescored exactly below.
            let _ = local_optimizer.minimize_with_model(
                &mut best_molecule,
                InterModel::exact(receptor),
                forcefield,
                bounds,
                &intra_pairs,
            );
        }

        // The search minimizes inter + intra, but Vina reports the affinity
        // from the intermolecular term alone: the ligand's internal energy
        // cancels against the unbound reference state. Rescore exactly so
        // grid interpolation error never reaches the reported numbers.
        let inter_energy = intermolecular_energy(&best_molecule, receptor, forcefield)?;
        let intra_energy = intramolecular_energy(&best_molecule, &intra_pairs, forcefield)?;

        let energy_components = vec![
            ("Intermolecular".to_string(), inter_energy),
            ("Intramolecular".to_string(), intra_energy),
        ];

        Ok(DockingResult {
            energy: inter_energy,
            molecule: best_molecule,
            energy_components,
            rmsd: None,
        })
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
    fn optimize(
        &self,
        ligand: &Molecule,
        receptor: &Molecule,
        forcefield: &dyn ForceField,
        center: Vector3<f64>,
        box_size: Vector3<f64>,
        max_steps: usize,
    ) -> Result<DockingResult, OptimizationError> {
        let grids = self.build_grids(ligand, receptor, forcefield, &center, &box_size)?;
        self.run(
            ligand,
            receptor,
            grids.as_ref(),
            forcefield,
            center,
            box_size,
            max_steps,
        )
    }

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
        // Receptor atoms too far from the box can never interact with any pose
        // inside it, so drop them once instead of testing them on every one of
        // the millions of energy evaluations to come.
        let pruned = prune_receptor(receptor, &center, &box_size, ligand);

        // One set of grids serves every run.
        let grids = self.build_grids(ligand, &pruned, forcefield, &center, &box_size)?;

        // Generate poses in parallel using rayon. If a master seed is set, each
        // run gets its own seed mixed from (master, run index), so the parallel
        // pose set stays reproducible and no two master seeds share runs.
        let master_seed = self.params.seed;
        let results: Result<Vec<DockingResult>, OptimizationError> = (0..num_poses)
            .into_par_iter()
            .map(|i| {
                let per_pose = MonteCarlo::with_params(MonteCarloParams {
                    seed: master_seed.map(|s| run_seed(s, i as u64)),
                    ..self.params.clone()
                });
                per_pose.run(
                    ligand,
                    &pruned,
                    grids.as_ref(),
                    forcefield,
                    center,
                    box_size,
                    max_steps,
                )
            })
            .collect();

        let mut results = results?;

        // Sort by energy
        results.sort_by(|a, b| a.energy.partial_cmp(&b.energy).unwrap());

        // Collapse near-identical poses. The runs are independent, so many of
        // them land in the same minimum; without this, "9 binding modes" can be
        // nine copies of one pose.
        let mut clustered: Vec<DockingResult> = Vec::new();
        for result in results {
            let duplicate = clustered
                .iter()
                .any(|kept| calculate_rmsd(&result.molecule, &kept.molecule) < MODE_RMSD_THRESHOLD);
            if !duplicate {
                clustered.push(result);
            }
        }
        let mut results = clustered;

        // Calculate RMSD to best pose
        if !results.is_empty() {
            let (first, rest) = results.split_at_mut(1);
            let best_molecule = &first[0].molecule;

            for result in rest {
                let rmsd = calculate_rmsd(&result.molecule, best_molecule);
                result.rmsd = Some(rmsd);
            }
        }

        Ok(results)
    }
}

/// Minimum RMSD between two reported binding modes, in Angstroms.
///
/// Vina uses the same idea: modes closer than this are the same answer found
/// twice, not two distinct hypotheses.
const MODE_RMSD_THRESHOLD: f64 = 1.0;

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

/// Sample a rotation axis uniformly over the sphere.
///
/// Sampling the components from a cube and normalizing would bias the axis
/// toward the eight corners; rejection sampling on the unit ball does not.
fn random_axis(rng: &mut StdRng) -> Unit<Vector3<f64>> {
    for _ in 0..32 {
        let v = Vector3::new(
            rng.gen::<f64>() * 2.0 - 1.0,
            rng.gen::<f64>() * 2.0 - 1.0,
            rng.gen::<f64>() * 2.0 - 1.0,
        );
        let norm_sq = v.norm_squared();
        if norm_sq > 1e-12 && norm_sq <= 1.0 {
            return Unit::new_normalize(v);
        }
    }
    Unit::new_unchecked(Vector3::new(0.0, 0.0, 1.0))
}

/// Sample a uniformly random orientation.
fn random_orientation(rng: &mut StdRng) -> UnitQuaternion<f64> {
    UnitQuaternion::from_axis_angle(&random_axis(rng), (rng.gen::<f64>() - 0.5) * 2.0 * PI)
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
