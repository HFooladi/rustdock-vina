//! AutoDock Vina's global search, ported from Vina 1.2.7.
//!
//! Each piece follows a Vina source file:
//!
//! | Here | Vina |
//! |---|---|
//! | [`Conf`], [`Model`] | `conf.h`, `model::set` |
//! | [`Searcher::bfgs`] | `bfgs.h`, `quasi_newton.cpp` |
//! | [`Searcher::mutate`] | `mutate.cpp` |
//! | [`Searcher::chain`] | `monte_carlo.cpp` |
//! | [`add_to_container`] | `coords.cpp` |
//! | [`Searcher::search`] | `parallel_mc.cpp`, `Vina::global_search` |
//!
//! Poses are evaluated on a coordinate-only model of the ligand and any
//! flexible side chains, so a trial pose costs a pass over a coordinate
//! array rather than a `Molecule` clone.

use std::collections::HashMap;
use std::f64::consts::PI;

use nalgebra::{DMatrix, DVector, Rotation3, Unit, UnitQuaternion, Vector3};
use rand::prelude::*;
use rand::rngs::StdRng;
use rayon::prelude::*;

use crate::atom::Atom;
use crate::forcefield::ForceField;
use crate::grid::{curl, GridSet};
use crate::molecule::{FlexResidues, Molecule};
use crate::optimization::OptimizationError;

/// Out-of-box penalty per Angstrom during the search (Vina's `slope`).
pub(crate) const SLOPE: f64 = 1e6;

/// The curl cap Vina uses outside the Monte Carlo moves (`authentic_v`):
/// large enough to leave realistic energies untouched.
pub(crate) const AUTHENTIC_CAP: f64 = 1000.0;

/// Curl caps for the two halves of the search energy (Vina's `v`): the
/// energy of every moving atom against the rigid receptor (`v[1]`), and the
/// pairs among moving atoms: the ligand's own, glue, ligand-flex and
/// flex-flex (`v[0]`, `v[2]`).
#[derive(Debug, Clone, Copy)]
pub(crate) struct Caps {
    pub(crate) receptor: f64,
    pub(crate) intra: f64,
}

/// Vina's `authentic_v`.
pub(crate) const AUTHENTIC: Caps = Caps {
    receptor: AUTHENTIC_CAP,
    intra: AUTHENTIC_CAP,
};

/// Weight of the linear attraction that closes an opened macrocycle
/// (Vina's `weight_glue`), and the distance it reaches to (`max_cutoff`).
const GLUE_WEIGHT: f64 = 50.0;
const GLUE_CUTOFF: f64 = 20.0;

/// RMSD within which two poses from different chains are merged.
const MERGE_RMSD: f64 = 2.0;

/// The search box as Vina grids it: a whole number of `spacing` cells
/// covering `size`, centred on `center`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct SearchBox {
    lo: Vector3<f64>,
    hi: Vector3<f64>,
}

impl SearchBox {
    pub(crate) fn new(center: &Vector3<f64>, size: &Vector3<f64>, spacing: f64) -> Self {
        let span = size.map(|s| ((s / spacing).ceil()).max(1.0) * spacing);
        let lo = center - span * 0.5;
        Self { lo, hi: lo + span }
    }

    /// `p` clamped into the box, the summed distance it lay outside, and the
    /// direction of the penalty gradient on each axis.
    fn clamp(&self, p: &Vector3<f64>) -> (Vector3<f64>, f64, Vector3<f64>) {
        let mut adjusted = *p;
        let mut miss = 0.0;
        let mut sign = Vector3::zeros();
        for i in 0..3 {
            if p[i] < self.lo[i] {
                adjusted[i] = self.lo[i];
                miss += self.lo[i] - p[i];
                sign[i] = -1.0;
            } else if p[i] > self.hi[i] {
                adjusted[i] = self.hi[i];
                miss += p[i] - self.hi[i];
                sign[i] = 1.0;
            }
        }
        (adjusted, miss, sign)
    }

    fn contains(&self, p: &Vector3<f64>, margin: f64) -> bool {
        (0..3).all(|i| p[i] >= self.lo[i] - margin && p[i] <= self.hi[i] + margin)
    }

    fn random_point(&self, rng: &mut StdRng) -> Vector3<f64> {
        Vector3::new(
            rng.gen_range(self.lo.x..=self.hi.x),
            rng.gen_range(self.lo.y..=self.hi.y),
            rng.gen_range(self.lo.z..=self.hi.z),
        )
    }
}

/// A pose: the ligand's rigid placement and torsion angles, and the flexible
/// side chains' torsion angles.
///
/// Ligand coordinates are `pivot + position + R(orientation) (p(torsions) -
/// pivot)`, where `p` applies the torsion angles to the input coordinates
/// and `pivot` is the input centroid. Side chains only turn about their
/// torsions. The identity [`Conf`] reproduces the input pose.
#[derive(Debug, Clone)]
pub(crate) struct Conf {
    position: Vector3<f64>,
    orientation: UnitQuaternion<f64>,
    torsions: Vec<f64>,
    flex: Vec<f64>,
}

impl Conf {
    fn identity(n_torsions: usize, n_flex: usize) -> Self {
        Self {
            position: Vector3::zeros(),
            orientation: UnitQuaternion::identity(),
            torsions: vec![0.0; n_torsions],
            flex: vec![0.0; n_flex],
        }
    }

    /// Step along `change` scaled by `factor` (Vina's `conf::increment`).
    /// The layout is Vina's: position, rotation, ligand torsions, then flex
    /// torsions. Rotation is a rotation vector applied on the left, so the
    /// orientation component of the gradient is the torque.
    fn increment(&mut self, change: &DVector<f64>, factor: f64) {
        self.position += factor * Vector3::new(change[0], change[1], change[2]);
        let rotation = factor * Vector3::new(change[3], change[4], change[5]);
        self.orientation = UnitQuaternion::from_scaled_axis(rotation) * self.orientation;
        let n = self.torsions.len();
        for (k, t) in self.torsions.iter_mut().enumerate() {
            *t = normalized_angle(*t + factor * change[6 + k]);
        }
        for (k, t) in self.flex.iter_mut().enumerate() {
            *t = normalized_angle(*t + factor * change[6 + n + k]);
        }
    }
}

struct TorsionAxis {
    /// Bond atoms: the axis runs from `a` to `b`.
    a: usize,
    b: usize,
    moving: Vec<usize>,
}

impl TorsionAxis {
    /// The axis of `molecule`'s torsion `t`, with atom indices passed
    /// through `index`; `None` if it names an atom that does not exist.
    fn of(molecule: &Molecule, t: &crate::molecule::Torsion, index: &[usize]) -> Option<Self> {
        let bond = molecule.bonds.get(t.bond_idx)?;
        let at = |i: usize| index.get(i).copied();
        Some(Self {
            a: at(bond.atom1_idx)?,
            b: at(bond.atom2_idx)?,
            moving: t
                .moving_atoms
                .iter()
                .map(|&i| at(i))
                .collect::<Option<_>>()?,
        })
    }

    /// Turn `out` by `angle` about this axis.
    fn apply(&self, angle: f64, out: &mut [Vector3<f64>]) {
        if angle == 0.0 {
            return;
        }
        let origin = out[self.a];
        let direction = out[self.b] - origin;
        if direction.norm_squared() < 1e-20 {
            return;
        }
        let rotation = Rotation3::from_axis_angle(&Unit::new_normalize(direction), angle);
        for &i in &self.moving {
            out[i] = rotation * (out[i] - origin) + origin;
        }
    }

    /// The derivative of the energy with respect to this torsion: the
    /// torque of the moving atoms' gradients about the axis.
    fn derivative(&self, coords: &[Vector3<f64>], grads: &[Vector3<f64>]) -> f64 {
        let origin = coords[self.a];
        let direction = coords[self.b] - origin;
        if direction.norm_squared() < 1e-20 {
            return 0.0;
        }
        let u = direction.normalize();
        self.moving
            .iter()
            .map(|&i| grads[i].dot(&u.cross(&(coords[i] - origin))))
            .sum()
    }
}

/// A coordinate-only model of what the search moves: the ligand and any
/// flexible side chains, enough to turn a [`Conf`] into coordinates and
/// score them.
///
/// Atoms are laid out as in Vina's `model`: the ligand's, then the side
/// chains' moving atoms, then their fixed atoms. Every atom before
/// `n_movable` is scored against the rigid receptor.
pub(crate) struct Model {
    pub(crate) atoms: Vec<Atom>,
    n_ligand: usize,
    n_movable: usize,
    base: Vec<Vector3<f64>>,
    pivot: Vector3<f64>,
    torsions: Vec<TorsionAxis>,
    flex_torsions: Vec<TorsionAxis>,
    /// The ligand's heavy atoms.
    heavy: Vec<usize>,
    /// Every moving heavy atom, the ligand's and the side chains'.
    movable_heavy: Vec<usize>,
    pairs: Vec<(usize, usize)>,
    /// Closure dummy - closure carbon pairs of an opened macrocycle.
    glue: Vec<(usize, usize)>,
    /// Ligand-flex and flex-flex pairs (Vina's `inter_pairs` and
    /// `other_pairs`).
    flex_pairs: Vec<(usize, usize)>,
    /// For each flex atom, from `n_ligand` on, its index in the flex
    /// residues.
    flex_source: Vec<usize>,
}

impl Model {
    pub(crate) fn new(ligand: &Molecule, flex: Option<&FlexResidues>) -> Self {
        let n_ligand = ligand.atoms.len();
        let ligand_index: Vec<usize> = (0..n_ligand).collect();
        let torsions = ligand
            .torsions
            .iter()
            .filter_map(|t| TorsionAxis::of(ligand, t, &ligand_index))
            .collect();
        let pivot = ligand.center().unwrap_or_else(|_| Vector3::zeros());
        let mut atoms = ligand.atoms.clone();

        // Flex atoms, moving ones first.
        let mut flex_source = Vec::new();
        let mut flex_torsions = Vec::new();
        let mut flex_pairs = Vec::new();
        if let Some(flex) = flex {
            flex_source.extend(flex.movable());
            flex_source.extend((0..flex.fixed.len()).filter(|&i| flex.fixed[i]));
            let mut index = vec![0; flex_source.len()];
            for (k, &i) in flex_source.iter().enumerate() {
                index[i] = n_ligand + k;
                atoms.push(flex.molecule.atoms[i].clone());
            }
            flex_torsions = flex
                .molecule
                .torsions
                .iter()
                .filter_map(|t| TorsionAxis::of(&flex.molecule, t, &index))
                .collect();
            flex_pairs = flex
                .internal_pairs()
                .into_iter()
                .map(|(i, j)| (index[i], index[j]))
                .collect();
            // Every ligand atom against every moving flex atom, but for
            // hydrogens and closure dummies, as Vina's `model::append`.
            let scored = |a: &Atom| !a.atom_type.is_hydrogen() && !a.atom_type.is_closure_dummy();
            for i in flex.movable() {
                for j in 0..n_ligand {
                    if scored(&flex.molecule.atoms[i]) && scored(&ligand.atoms[j]) {
                        flex_pairs.push((j, index[i]));
                    }
                }
            }
        }
        let n_movable = n_ligand + flex.map_or(0, |f| f.movable().count());
        let heavy = (0..n_ligand)
            .filter(|&i| !atoms[i].atom_type.is_hydrogen())
            .collect();
        let movable_heavy = (0..n_movable)
            .filter(|&i| !atoms[i].atom_type.is_hydrogen())
            .collect();
        Self {
            base: atoms.iter().map(|a| a.coordinates).collect(),
            atoms,
            n_ligand,
            n_movable,
            pivot,
            torsions,
            flex_torsions,
            heavy,
            movable_heavy,
            pairs: ligand.intramolecular_pairs(),
            glue: ligand.glue_pairs(),
            flex_pairs,
            flex_source,
        }
    }

    pub(crate) fn n_torsions(&self) -> usize {
        self.torsions.len()
    }

    pub(crate) fn n_flex_torsions(&self) -> usize {
        self.flex_torsions.len()
    }

    /// Number of atoms the search moves, Vina's `num_movable_atoms`.
    pub(crate) fn n_movable(&self) -> usize {
        self.n_movable
    }

    fn identity(&self) -> Conf {
        Conf::identity(self.n_torsions(), self.n_flex_torsions())
    }

    /// Write the coordinates of `conf` into `out`.
    fn coords(&self, conf: &Conf, out: &mut [Vector3<f64>]) {
        out.copy_from_slice(&self.base);
        // Torsion moving sets are nested or disjoint, so applying each about
        // its current axis gives the same pose in any order.
        for (axis, &angle) in self.torsions.iter().zip(&conf.torsions) {
            axis.apply(angle, out);
        }
        let rotation = conf.orientation.to_rotation_matrix();
        let center = self.pivot + conf.position;
        for p in &mut out[..self.n_ligand] {
            *p = rotation * (*p - self.pivot) + center;
        }
        for (axis, &angle) in self.flex_torsions.iter().zip(&conf.flex) {
            axis.apply(angle, out);
        }
    }

    /// Chain rule from per-atom gradients at `coords` (the coordinates of
    /// `conf`) to the [`Conf`] tangent space: force and torque on the ligand
    /// about its rotation centre, then one torque projection per torsion.
    fn conf_gradient(
        &self,
        conf: &Conf,
        coords: &[Vector3<f64>],
        grads: &[Vector3<f64>],
    ) -> DVector<f64> {
        let n = self.torsions.len();
        let mut g = DVector::zeros(6 + n + self.flex_torsions.len());
        let center = self.pivot + conf.position;
        let mut force = Vector3::zeros();
        let mut torque = Vector3::zeros();
        for (p, d) in coords.iter().zip(grads).take(self.n_ligand) {
            force += d;
            torque += (p - center).cross(d);
        }
        for k in 0..3 {
            g[k] = force[k];
            g[3 + k] = torque[k];
        }
        for (k, axis) in self.torsions.iter().enumerate() {
            g[6 + k] = axis.derivative(coords, grads);
        }
        for (k, axis) in self.flex_torsions.iter().enumerate() {
            g[6 + n + k] = axis.derivative(coords, grads);
        }
        g
    }

    /// Root-mean-square distance of the ligand's heavy atoms from `center`
    /// (Vina's `gyration_radius`).
    fn gyration_radius(&self, coords: &[Vector3<f64>], center: &Vector3<f64>) -> f64 {
        if self.heavy.is_empty() {
            return 0.0;
        }
        let sum: f64 = self
            .heavy
            .iter()
            .map(|&i| (coords[i] - center).norm_squared())
            .sum();
        (sum / self.heavy.len() as f64).sqrt()
    }

    /// The ligand with its atoms moved to `coords`.
    pub(crate) fn to_molecule(
        &self,
        template: &Molecule,
        conf: &Conf,
        coords: &[Vector3<f64>],
    ) -> Molecule {
        let mut molecule = template.clone();
        for (atom, p) in molecule.atoms.iter_mut().zip(coords) {
            atom.coordinates = *p;
        }
        for (torsion, angle) in molecule.torsions.iter_mut().zip(&conf.torsions) {
            torsion.angle = normalized_angle(torsion.angle + angle);
        }
        molecule
    }

    /// The flexible side chains with their atoms moved to `coords`.
    pub(crate) fn to_flex(
        &self,
        template: &FlexResidues,
        conf: &Conf,
        coords: &[Vector3<f64>],
    ) -> FlexResidues {
        let mut flex = template.clone();
        for (&i, p) in self.flex_source.iter().zip(&coords[self.n_ligand..]) {
            flex.molecule.atoms[i].coordinates = *p;
        }
        for (torsion, angle) in flex.molecule.torsions.iter_mut().zip(&conf.flex) {
            torsion.angle = normalized_angle(torsion.angle + angle);
        }
        flex
    }
}

/// Intramolecular pair energies tabulated against squared distance, one
/// table per pair of ligand atom classes (Vina's `precalculate`).
///
/// Only valid for force fields whose pair energy depends on nothing but the
/// atoms' types, `xs` flags and distance (those with a
/// [`grid_cutoff`](ForceField::grid_cutoff)).
struct PairTables {
    tables: Vec<Vec<f64>>,
    /// Table index for each pair it was built for.
    pair_table: Vec<usize>,
    cutoff2: f64,
}

/// Samples per square Angstrom (Vina uses 32).
const PAIR_TABLE_FACTOR: f64 = 64.0;

impl PairTables {
    fn build(
        atoms: &[Atom],
        pairs: &[(usize, usize)],
        forcefield: &dyn ForceField,
        cutoff: f64,
    ) -> Result<Self, OptimizationError> {
        let key = |a: &Atom| (a.atom_type, a.is_hydrophobic, a.is_hbond_donor);
        let cutoff2 = cutoff * cutoff;
        let samples = (cutoff2 * PAIR_TABLE_FACTOR) as usize + 2;
        let mut index: HashMap<_, usize> = HashMap::new();
        let mut tables = Vec::new();
        let mut pair_table = Vec::with_capacity(pairs.len());
        for &(i, j) in pairs {
            let (a, b) = (&atoms[i], &atoms[j]);
            let (ka, kb) = (key(a), key(b));
            // Pair energy is symmetric, so either order finds the table.
            let found = index.get(&(ka, kb)).or_else(|| index.get(&(kb, ka)));
            let t = match found {
                Some(&t) => t,
                None => {
                    let mut table = Vec::with_capacity(samples);
                    for n in 0..samples {
                        let r = (n as f64 / PAIR_TABLE_FACTOR).sqrt().clamp(0.01, cutoff);
                        table.push(forcefield.atom_pair_energy(a, b, r)?);
                    }
                    tables.push(table);
                    index.insert((ka, kb), tables.len() - 1);
                    tables.len() - 1
                }
            };
            pair_table.push(t);
        }
        Ok(Self {
            tables,
            pair_table,
            cutoff2,
        })
    }

    /// Energy of pair number `pair` at squared distance `r2`, and its
    /// derivative with respect to `r2` (exact for the interpolant).
    #[inline]
    fn eval(&self, pair: usize, r2: f64) -> (f64, f64) {
        let table = &self.tables[self.pair_table[pair]];
        let x = r2 * PAIR_TABLE_FACTOR;
        let k = (x as usize).min(table.len() - 2);
        let t = x - k as f64;
        let (e0, e1) = (table[k], table[k + 1]);
        (e0 + t * (e1 - e0), (e1 - e0) * PAIR_TABLE_FACTOR)
    }
}

/// Energy against the rigid receptor as the search sees it.
pub(crate) trait InterField: Sync {
    /// Energy of the moving atoms at `coords`, adding each atom's gradient
    /// into `grads`. Atoms outside the box cost `slope` per Angstrom; `cap`
    /// softens each atom's positive energy.
    fn eval(
        &self,
        atoms: &[Atom],
        coords: &[Vector3<f64>],
        slope: f64,
        cap: f64,
        grads: &mut [Vector3<f64>],
    ) -> Result<f64, OptimizationError>;
}

/// Precomputed affinity grids covering the search box (Vina's `cache`).
pub(crate) struct GridField<'a>(pub(crate) &'a GridSet);

impl InterField for GridField<'_> {
    fn eval(
        &self,
        _atoms: &[Atom],
        coords: &[Vector3<f64>],
        slope: f64,
        cap: f64,
        grads: &mut [Vector3<f64>],
    ) -> Result<f64, OptimizationError> {
        Ok(self.0.search_energy(coords, slope, cap, grads))
    }
}

/// Exact pairwise energy against the receptor (Vina's `non_cache`): each
/// moving atom is clamped into the box, summed over receptor atoms within
/// the cutoff, curled, then charged the out-of-box penalty.
pub(crate) struct ExactField<'a> {
    receptor: Vec<&'a Atom>,
    forcefield: &'a dyn ForceField,
    search_box: SearchBox,
    cutoff: f64,
    cells: HashMap<(i64, i64, i64), Vec<usize>>,
}

impl<'a> ExactField<'a> {
    pub(crate) fn new(
        receptor: &'a Molecule,
        forcefield: &'a dyn ForceField,
        search_box: SearchBox,
    ) -> Self {
        let cutoff = forcefield
            .grid_cutoff()
            .unwrap_or(crate::forcefield::vina::CUTOFF);
        let receptor: Vec<&Atom> = receptor
            .atoms
            .iter()
            .filter(|a| !a.atom_type.is_hydrogen())
            .collect();
        let mut cells: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
        for (i, atom) in receptor.iter().enumerate() {
            cells
                .entry(cell_of(&atom.coordinates, cutoff))
                .or_default()
                .push(i);
        }
        Self {
            receptor,
            forcefield,
            search_box,
            cutoff,
            cells,
        }
    }
}

fn cell_of(p: &Vector3<f64>, size: f64) -> (i64, i64, i64) {
    (
        (p.x / size).floor() as i64,
        (p.y / size).floor() as i64,
        (p.z / size).floor() as i64,
    )
}

impl InterField for ExactField<'_> {
    fn eval(
        &self,
        atoms: &[Atom],
        coords: &[Vector3<f64>],
        slope: f64,
        cap: f64,
        grads: &mut [Vector3<f64>],
    ) -> Result<f64, OptimizationError> {
        let mut total = 0.0;
        let cutoff2 = self.cutoff * self.cutoff;
        for ((atom, p), grad) in atoms.iter().zip(coords).zip(grads) {
            if atom.atom_type.is_hydrogen() {
                continue;
            }
            let (adjusted, miss, sign) = self.search_box.clamp(p);
            let mut probe = atom.clone();
            probe.coordinates = adjusted;
            let mut e = 0.0;
            let mut d = Vector3::zeros();
            let (cx, cy, cz) = cell_of(&adjusted, self.cutoff);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let Some(members) = self.cells.get(&(cx + dx, cy + dy, cz + dz)) else {
                            continue;
                        };
                        for &j in members {
                            let rec = self.receptor[j];
                            let delta = adjusted - rec.coordinates;
                            let r2 = delta.norm_squared();
                            if r2 >= cutoff2 || r2 < 1e-4 {
                                continue;
                            }
                            let r = r2.sqrt();
                            e += self.forcefield.atom_pair_energy(&probe, rec, r)?;
                            d +=
                                delta * (self.forcefield.atom_pair_derivative(&probe, rec, r)? / r);
                        }
                    }
                }
            }
            curl(&mut e, &mut d, cap);
            total += e + slope * miss;
            *grad += d + slope * sign;
        }
        Ok(total)
    }
}

/// Search settings, as Vina sets them in `Vina::global_search`.
#[derive(Debug, Clone)]
pub(crate) struct SearchSettings {
    pub(crate) chains: usize,
    pub(crate) global_steps: usize,
    pub(crate) local_steps: usize,
    pub(crate) num_saved: usize,
    pub(crate) min_rmsd: f64,
    pub(crate) temperature: f64,
    pub(crate) mutation_amplitude: f64,
    pub(crate) hunt_caps: Caps,
    pub(crate) randomize: bool,
}

/// Vina's Monte Carlo step count per chain:
/// `70 * 3 * (50 + atoms + 10 * dof) / 2`, over the moving atoms and the
/// degrees of freedom (6 + torsions, the side chains' included).
pub(crate) fn vina_global_steps(movable_atoms: usize, torsions: usize) -> usize {
    let dof = 6 + torsions;
    let heuristic = movable_atoms + 10 * dof;
    70 * 3 * (50 + heuristic) / 2
}

/// Vina's BFGS step count: `(25 + moving atoms) / 3`.
pub(crate) fn vina_local_steps(movable_atoms: usize) -> usize {
    (25 + movable_atoms) / 3
}

/// A pose found by the search.
#[derive(Debug, Clone)]
pub(crate) struct Pose {
    pub(crate) conf: Conf,
    pub(crate) e: f64,
    heavy: Vec<Vector3<f64>>,
}

/// Reusable buffers for evaluating poses.
struct Scratch {
    coords: Vec<Vector3<f64>>,
    grads: Vec<Vector3<f64>>,
}

/// The search over one ligand in one box.
pub(crate) struct Searcher<'a> {
    pub(crate) model: Model,
    forcefield: &'a dyn ForceField,
    search_box: SearchBox,
    pair_tables: Option<PairTables>,
    flex_tables: Option<PairTables>,
}

impl<'a> Searcher<'a> {
    pub(crate) fn new(
        ligand: &Molecule,
        flex: Option<&FlexResidues>,
        forcefield: &'a dyn ForceField,
        search_box: SearchBox,
    ) -> Result<Self, OptimizationError> {
        let model = Model::new(ligand, flex);
        let tables = |pairs: &[(usize, usize)]| match forcefield.grid_cutoff() {
            Some(cutoff) => PairTables::build(&model.atoms, pairs, forcefield, cutoff).map(Some),
            None => Ok(None),
        };
        let pair_tables = tables(&model.pairs)?;
        let flex_tables = tables(&model.flex_pairs)?;
        Ok(Self {
            model,
            forcefield,
            search_box,
            pair_tables,
            flex_tables,
        })
    }

    fn scratch(&self) -> Scratch {
        let n = self.model.base.len();
        Scratch {
            coords: vec![Vector3::zeros(); n],
            grads: vec![Vector3::zeros(); n],
        }
    }

    /// Search energy (inter + intra) of `conf` and its [`Conf`] gradient
    /// (Vina's `model::eval_deriv`): every moving atom against the rigid
    /// receptor, then the pairs among them.
    fn eval_deriv(
        &self,
        field: &dyn InterField,
        conf: &Conf,
        slope: f64,
        caps: Caps,
        s: &mut Scratch,
    ) -> Result<(f64, DVector<f64>), OptimizationError> {
        self.model.coords(conf, &mut s.coords);
        s.grads.iter_mut().for_each(|g| *g = Vector3::zeros());
        let n = self.model.n_movable;
        let mut e = field.eval(
            &self.model.atoms[..n],
            &s.coords[..n],
            slope,
            caps.receptor,
            &mut s.grads[..n],
        )?;
        e += self.intra(&s.coords, caps.intra, Some(&mut s.grads))?;
        e += self.flex_pairs(&s.coords, caps.intra, Some(&mut s.grads))?;
        e += self.glue(&s.coords, caps.intra, &mut s.grads)?;
        Ok((e, self.model.conf_gradient(conf, &s.coords, &s.grads)))
    }

    /// The pull that closes an opened macrocycle: each glue pair's pair
    /// energy plus [`GLUE_WEIGHT`] per Angstrom of separation, curled with
    /// `cap`. Part of the search energy only, never of a reported one.
    fn glue(
        &self,
        coords: &[Vector3<f64>],
        cap: f64,
        grads: &mut [Vector3<f64>],
    ) -> Result<f64, OptimizationError> {
        let atoms = &self.model.atoms;
        let mut total = 0.0;
        for &(i, j) in &self.model.glue {
            let delta = coords[i] - coords[j];
            let r = delta.norm();
            if r >= GLUE_CUTOFF || r <= 1e-6 {
                continue;
            }
            let mut e =
                self.forcefield.atom_pair_energy(&atoms[i], &atoms[j], r)? + GLUE_WEIGHT * r;
            let de_dr = self
                .forcefield
                .atom_pair_derivative(&atoms[i], &atoms[j], r)?
                + GLUE_WEIGHT;
            let mut d = delta * (de_dr / r);
            curl(&mut e, &mut d, cap);
            total += e;
            grads[i] += d;
            grads[j] -= d;
        }
        Ok(total)
    }

    /// The ligand's internal energy over its interacting pairs, each pair's
    /// positive energy curled with `cap`.
    fn intra(
        &self,
        coords: &[Vector3<f64>],
        cap: f64,
        grads: Option<&mut [Vector3<f64>]>,
    ) -> Result<f64, OptimizationError> {
        let tables = self.pair_tables.as_ref();
        self.pair_energy(&self.model.pairs, tables, coords, cap, grads)
    }

    /// The flexible side chains' energy with the ligand and among
    /// themselves, each pair curled with `cap`.
    fn flex_pairs(
        &self,
        coords: &[Vector3<f64>],
        cap: f64,
        grads: Option<&mut [Vector3<f64>]>,
    ) -> Result<f64, OptimizationError> {
        let tables = self.flex_tables.as_ref();
        self.pair_energy(&self.model.flex_pairs, tables, coords, cap, grads)
    }

    /// The energy of `pairs`, from `tables` (built for them) when there
    /// are, each pair's positive energy curled with `cap`.
    fn pair_energy(
        &self,
        pairs: &[(usize, usize)],
        tables: Option<&PairTables>,
        coords: &[Vector3<f64>],
        cap: f64,
        mut grads: Option<&mut [Vector3<f64>]>,
    ) -> Result<f64, OptimizationError> {
        let atoms = &self.model.atoms;
        let mut total = 0.0;
        if let Some(tables) = tables {
            for (n, &(i, j)) in pairs.iter().enumerate() {
                let delta = coords[i] - coords[j];
                let r2 = delta.norm_squared();
                if r2 >= tables.cutoff2 || r2 <= 1e-4 {
                    continue;
                }
                let (mut e, de_dr2) = tables.eval(n, r2);
                let mut d = delta * (2.0 * de_dr2);
                curl(&mut e, &mut d, cap);
                total += e;
                if let Some(g) = grads.as_deref_mut() {
                    g[i] += d;
                    g[j] -= d;
                }
            }
            return Ok(total);
        }
        for &(i, j) in pairs {
            let delta = coords[i] - coords[j];
            let r = delta.norm();
            if r <= 0.01 {
                continue;
            }
            let mut e = self.forcefield.atom_pair_energy(&atoms[i], &atoms[j], r)?;
            if e == 0.0 {
                continue;
            }
            let mut d = delta
                * (self
                    .forcefield
                    .atom_pair_derivative(&atoms[i], &atoms[j], r)?
                    / r);
            curl(&mut e, &mut d, cap);
            total += e;
            if let Some(g) = grads.as_deref_mut() {
                g[i] += d;
                g[j] -= d;
            }
        }
        Ok(total)
    }

    /// Vina's BFGS (`bfgs.h`), minimizing in place; returns the final energy.
    ///
    /// Kept faithful down to its quirks: the convergence test reads the
    /// previous gradient, a failed line search still takes its last trial,
    /// and the result reverts to the start if it ended higher.
    fn bfgs(
        &self,
        field: &dyn InterField,
        x: &mut Conf,
        slope: f64,
        caps: Caps,
        max_steps: usize,
        s: &mut Scratch,
    ) -> Result<f64, OptimizationError> {
        let n = 6 + self.model.n_torsions() + self.model.n_flex_torsions();
        let mut h = DMatrix::<f64>::identity(n, n);
        let (mut f0, mut g) = self.eval_deriv(field, x, slope, caps, s)?;
        let (f_orig, x_orig) = (f0, x.clone());

        for step in 0..max_steps {
            let p = -(&h * &g);
            let pg = p.dot(&g);

            // Backtracking line search: alpha from 1, halved up to 10 times.
            let mut alpha = 1.0;
            let mut x_new = x.clone();
            let mut f1 = f0;
            let mut g_new = g.clone();
            for _ in 0..10 {
                x_new = x.clone();
                x_new.increment(&p, alpha);
                (f1, g_new) = self.eval_deriv(field, &x_new, slope, caps, s)?;
                if f1 - f0 < 1e-4 * alpha * pg {
                    break;
                }
                alpha *= 0.5;
            }

            let y = &g_new - &g;
            f0 = f1;
            *x = x_new;
            // Vina's test, which also stops on NaN.
            let norm = g.norm();
            if norm < 1e-5 || norm.is_nan() {
                break;
            }
            g = g_new;

            if step == 0 {
                let yy = y.dot(&y);
                if yy.abs() > f64::EPSILON {
                    h = DMatrix::identity(n, n) * (alpha * y.dot(&p) / yy);
                }
            }
            bfgs_update(&mut h, &p, &y, alpha);
        }
        if f0 > f_orig || f0.is_nan() {
            f0 = f_orig;
            *x = x_orig;
        }
        Ok(f0)
    }

    /// Change exactly one thing (Vina's `mutate_conf`): the position by up to
    /// `amplitude` Angstroms, the orientation by up to `amplitude` over the
    /// gyration radius, or one ligand or side-chain torsion to a uniformly
    /// random angle.
    fn mutate(&self, conf: &mut Conf, amplitude: f64, rng: &mut StdRng, s: &mut Scratch) {
        let n = self.model.n_torsions();
        let entities = 2 + n + self.model.n_flex_torsions();
        let which = rng.gen_range(0..entities);
        match which {
            0 => conf.position += amplitude * random_inside_sphere(rng),
            1 => {
                self.model.coords(conf, &mut s.coords);
                let center = self.model.pivot + conf.position;
                let gr = self.model.gyration_radius(&s.coords, &center);
                if gr > f64::EPSILON {
                    let rotation = (amplitude / gr) * random_inside_sphere(rng);
                    conf.orientation =
                        UnitQuaternion::from_scaled_axis(rotation) * conf.orientation;
                }
            }
            k if k < 2 + n => conf.torsions[k - 2] = rng.gen_range(-PI..PI),
            k => conf.flex[k - 2 - n] = rng.gen_range(-PI..PI),
        }
    }

    /// A uniformly random pose: position in the box, orientation and every
    /// torsion, the side chains' included, uniformly random (Vina's
    /// `conf::randomize`).
    fn random_conf(&self, rng: &mut StdRng) -> Conf {
        let target = self.search_box.random_point(rng);
        let orientation = UnitQuaternion::from_quaternion(nalgebra::Quaternion::new(
            standard_normal(rng),
            standard_normal(rng),
            standard_normal(rng),
            standard_normal(rng),
        ));
        Conf {
            position: target - self.model.pivot,
            orientation,
            torsions: (0..self.model.n_torsions())
                .map(|_| rng.gen_range(-PI..PI))
                .collect(),
            flex: (0..self.model.n_flex_torsions())
                .map(|_| rng.gen_range(-PI..PI))
                .collect(),
        }
    }

    /// The moving heavy atoms' coordinates, which poses are clustered on
    /// (Vina's `get_heavy_atom_movable_coords`).
    fn heavy_coords(&self, conf: &Conf, s: &mut Scratch) -> Vec<Vector3<f64>> {
        self.model.coords(conf, &mut s.coords);
        self.model
            .movable_heavy
            .iter()
            .map(|&i| s.coords[i])
            .collect()
    }

    /// One Monte Carlo chain (Vina's `monte_carlo::operator()`): mutate,
    /// minimize, Metropolis, and keep the distinct low minima it visits.
    fn chain(
        &self,
        field: &dyn InterField,
        settings: &SearchSettings,
        seed: u64,
    ) -> Result<Vec<Pose>, OptimizationError> {
        let mut rng = StdRng::seed_from_u64(seed);
        let mut s = self.scratch();
        let mut current = if settings.randomize {
            self.random_conf(&mut rng)
        } else {
            self.model.identity()
        };
        let mut current_e = 0.0;
        let mut best_e = f64::INFINITY;
        let mut out: Vec<Pose> = Vec::new();

        for step in 0..settings.global_steps {
            let mut candidate = current.clone();
            self.mutate(
                &mut candidate,
                settings.mutation_amplitude,
                &mut rng,
                &mut s,
            );
            let candidate_e = self.bfgs(
                field,
                &mut candidate,
                SLOPE,
                settings.hunt_caps,
                settings.local_steps,
                &mut s,
            )?;
            if step == 0
                || metropolis_accept(current_e, candidate_e, settings.temperature, &mut rng)
            {
                current = candidate;
                current_e = candidate_e;
                if current_e < best_e || out.len() < settings.num_saved {
                    current_e = self.bfgs(
                        field,
                        &mut current,
                        SLOPE,
                        AUTHENTIC,
                        settings.local_steps,
                        &mut s,
                    )?;
                    let pose = Pose {
                        heavy: self.heavy_coords(&current, &mut s),
                        conf: current.clone(),
                        e: current_e,
                    };
                    add_to_container(&mut out, pose, settings.min_rmsd, settings.num_saved);
                    best_e = best_e.min(current_e);
                }
            }
        }
        Ok(out)
    }

    /// Vina's global search: independent chains in parallel, merged, then
    /// each pose refined and rescored on exact energies.
    ///
    /// Returns poses best first, with their search energies (inter + intra).
    pub(crate) fn search(
        &self,
        field: &dyn InterField,
        exact: &dyn InterField,
        settings: &SearchSettings,
        chain_seed: impl Fn(usize) -> u64 + Sync,
    ) -> Result<Vec<Pose>, OptimizationError> {
        let chains: Vec<Vec<Pose>> = (0..settings.chains)
            .into_par_iter()
            .map(|i| self.chain(field, settings, chain_seed(i)))
            .collect::<Result<_, _>>()?;

        let mut merged: Vec<Pose> = Vec::new();
        for pose in chains.into_iter().flatten() {
            add_to_container(&mut merged, pose, MERGE_RMSD, settings.num_saved);
        }
        let mut poses: Vec<Pose> = Vec::new();
        let n = merged.len();
        for pose in merged {
            add_to_container(&mut poses, pose, settings.min_rmsd, n);
        }

        // Refine on exact energies, raising the out-of-box slope until the
        // pose is inside, then rescore with the search slope.
        let mut s = self.scratch();
        for pose in &mut poses {
            for p in 0..5 {
                let slope = 100.0 * 10f64.powi(2 * p);
                self.bfgs(
                    exact,
                    &mut pose.conf,
                    slope,
                    AUTHENTIC,
                    settings.local_steps,
                    &mut s,
                )?;
                if self.within_box(&pose.conf, &mut s) {
                    break;
                }
            }
            pose.e = self
                .eval_deriv(exact, &pose.conf, SLOPE, AUTHENTIC, &mut s)?
                .0;
            pose.heavy = self.heavy_coords(&pose.conf, &mut s);
        }
        poses.sort_by(|a, b| a.e.total_cmp(&b.e));
        Ok(poses)
    }

    fn within_box(&self, conf: &Conf, s: &mut Scratch) -> bool {
        self.model.coords(conf, &mut s.coords);
        self.model
            .movable_heavy
            .iter()
            .all(|&i| self.search_box.contains(&s.coords[i], 1e-4))
    }

    /// The atoms scored against the rigid receptor, in the order the
    /// field sees them.
    pub(crate) fn movable_atoms(&self) -> &[Atom] {
        &self.model.atoms[..self.model.n_movable]
    }

    /// Coordinates of `conf`.
    pub(crate) fn coords_of(&self, conf: &Conf) -> Vec<Vector3<f64>> {
        let mut out = vec![Vector3::zeros(); self.model.base.len()];
        self.model.coords(conf, &mut out);
        out
    }

    /// The ligand's internal energy at `coords`, uncapped as Vina reports it.
    pub(crate) fn intra_energy(&self, coords: &[Vector3<f64>]) -> Result<f64, OptimizationError> {
        self.intra(coords, AUTHENTIC_CAP, None)
    }
}

/// Vina's BFGS inverse-Hessian update with step `alpha * p` and gradient
/// change `y`; skipped when the curvature condition fails.
fn bfgs_update(h: &mut DMatrix<f64>, p: &DVector<f64>, y: &DVector<f64>, alpha: f64) {
    let yp = y.dot(p);
    if alpha * yp < f64::EPSILON {
        return;
    }
    let minus_hy = -(&*h * y);
    let yhy = -y.dot(&minus_hy);
    let r = 1.0 / (alpha * yp);
    let n = p.len();
    for i in 0..n {
        for j in i..n {
            let v = h[(i, j)]
                + alpha * r * (minus_hy[i] * p[j] + minus_hy[j] * p[i])
                + alpha * alpha * (r * r * yhy + r) * p[i] * p[j];
            h[(i, j)] = v;
            h[(j, i)] = v;
        }
    }
}

/// Add `pose` to an energy-sorted container (Vina's `add_to_output_container`):
/// a pose within `min_rmsd` of a kept one replaces it if lower; otherwise it
/// is added, or replaces the worst when the container is full and it is
/// better.
fn add_to_container(out: &mut Vec<Pose>, pose: Pose, min_rmsd: f64, max_size: usize) {
    let closest = out
        .iter()
        .enumerate()
        .map(|(i, kept)| (i, rmsd(&pose.heavy, &kept.heavy)))
        .min_by(|a, b| a.1.total_cmp(&b.1));
    match closest {
        Some((i, d)) if d < min_rmsd => {
            if pose.e < out[i].e {
                out[i] = pose;
            }
        }
        _ => {
            if out.len() < max_size {
                out.push(pose);
            } else if out.last().is_some_and(|worst| pose.e < worst.e) {
                *out.last_mut().unwrap() = pose;
            }
        }
    }
    out.sort_by(|a, b| a.e.total_cmp(&b.e));
}

fn rmsd(a: &[Vector3<f64>], b: &[Vector3<f64>]) -> f64 {
    if a.is_empty() || a.len() != b.len() {
        return f64::INFINITY;
    }
    let sum: f64 = a.iter().zip(b).map(|(p, q)| (p - q).norm_squared()).sum();
    (sum / a.len() as f64).sqrt()
}

fn metropolis_accept(old: f64, new: f64, temperature: f64, rng: &mut StdRng) -> bool {
    new < old || rng.gen::<f64>() < ((old - new) / temperature).exp()
}

/// A point uniformly distributed inside the unit sphere.
fn random_inside_sphere(rng: &mut StdRng) -> Vector3<f64> {
    loop {
        let v = Vector3::new(
            rng.gen_range(-1.0..1.0),
            rng.gen_range(-1.0..1.0),
            rng.gen_range(-1.0..1.0),
        );
        if v.norm_squared() < 1.0 {
            return v;
        }
    }
}

/// A standard normal sample (Box-Muller).
fn standard_normal(rng: &mut StdRng) -> f64 {
    let u1: f64 = rng.gen_range(f64::EPSILON..1.0);
    let u2: f64 = rng.gen();
    (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
}

/// Wrap an angle into `(-pi, pi]`.
fn normalized_angle(a: f64) -> f64 {
    let wrapped = (a + PI).rem_euclid(2.0 * PI) - PI;
    if wrapped == -PI {
        PI
    } else {
        wrapped
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::vina::VinaForceField;
    use crate::io::parse_pdbqt;
    use std::path::PathBuf;

    fn data(name: &str) -> Molecule {
        parse_pdbqt(
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("benchmark_data")
                .join(name),
        )
        .unwrap()
    }

    /// The analytic [`Conf`] gradient must match finite differences along
    /// each tangent direction, for both fields, on a flexible ligand.
    #[test]
    fn test_conf_gradient_matches_finite_difference() {
        let receptor = data("receptor.pdbqt");
        let ligand = data("ligand.pdbqt");
        let center = ligand.center().unwrap();
        check_conf_gradient(&receptor, &ligand, None, center, Vector3::repeat(20.0), 6);
    }

    /// The same with three flexible side chains, whose torsions follow the
    /// ligand's in the tangent space.
    #[test]
    fn test_flex_conf_gradient_matches_finite_difference() {
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/parity_data/5SAK_ZRY");
        let receptor = parse_pdbqt(dir.join("flex_rigid.pdbqt")).unwrap();
        let ligand = parse_pdbqt(dir.join("ligand_xtal.pdbqt")).unwrap();
        let flex = crate::io::parse_flex_pdbqt(dir.join("flex.pdbqt"), &receptor).unwrap();
        let center = ligand.center().unwrap();
        let torsions = ligand.torsions.len() + flex.molecule.torsions.len();
        check_conf_gradient(
            &receptor,
            &ligand,
            Some(&flex),
            center,
            Vector3::repeat(24.0),
            torsions,
        );
    }

    /// Vina's energy steps at the 8 A cutoff, so a difference whose step
    /// straddles a pair crossing it is meaningless; those components are
    /// recognised by disagreeing one-sided differences and skipped.
    fn check_conf_gradient(
        receptor: &Molecule,
        ligand: &Molecule,
        flex: Option<&FlexResidues>,
        center: Vector3<f64>,
        size: Vector3<f64>,
        torsions: usize,
    ) {
        let ff = VinaForceField::default();
        let search_box = SearchBox::new(&center, &size, 0.375);
        let searcher = Searcher::new(ligand, flex, &ff, search_box).unwrap();
        let mut movable = Molecule::new("movable");
        movable.atoms = searcher.movable_atoms().to_vec();
        let grids = GridSet::build(&movable, receptor, &ff, &center, &size, 0.0, 0.375)
            .unwrap()
            .unwrap();
        let exact = ExactField::new(receptor, &ff, search_box);
        let mut rng = StdRng::seed_from_u64(3);
        let mut s = searcher.scratch();
        let n_torsions = searcher.model.n_torsions();
        let n = 6 + n_torsions + searcher.model.n_flex_torsions();
        assert_eq!(n, 6 + torsions);

        let (mut checked, mut skipped) = (0, 0);
        let fields: [&dyn InterField; 2] = [&GridField(&grids), &exact];
        for field in fields {
            for _ in 0..4 {
                let mut conf = searcher.model.identity();
                for t in conf.torsions.iter_mut().chain(conf.flex.iter_mut()) {
                    *t = rng.gen_range(-0.3..0.3);
                }
                conf.position = 0.3 * random_inside_sphere(&mut rng);
                let (f0, g) = searcher
                    .eval_deriv(field, &conf, SLOPE, AUTHENTIC, &mut s)
                    .unwrap();
                for k in 0..n {
                    let mut e = DVector::zeros(n);
                    e[k] = 1.0;
                    let h = 1e-6;
                    let mut plus = conf.clone();
                    plus.increment(&e, h);
                    let mut minus = conf.clone();
                    minus.increment(&e, -h);
                    let mut eval = |c: &Conf| {
                        searcher
                            .eval_deriv(field, c, SLOPE, AUTHENTIC, &mut s)
                            .unwrap()
                            .0
                    };
                    let (fp, fm) = (eval(&plus), eval(&minus));
                    let (forward, backward) = ((fp - f0) / h, (f0 - fm) / h);
                    if (forward - backward).abs() > 1e-2 * (1.0 + forward.abs()) {
                        skipped += 1;
                        continue;
                    }
                    checked += 1;
                    let numeric = (fp - fm) / (2.0 * h);
                    assert!(
                        (numeric - g[k]).abs() < 1e-3 * (1.0 + numeric.abs()),
                        "component {k}: analytic {} vs numeric {numeric}",
                        g[k]
                    );
                }
            }
        }
        assert!(
            skipped * 10 < checked,
            "{skipped} components skipped, {checked} checked"
        );
    }

    /// Tabulated intramolecular energy matches direct evaluation.
    #[test]
    fn test_pair_tables_match_direct_intra() {
        let ligand = data("ligand.pdbqt");
        let ff = VinaForceField::default();
        let search_box = SearchBox::new(&Vector3::zeros(), &Vector3::repeat(20.0), 0.375);
        let tabled = Searcher::new(&ligand, None, &ff, search_box).unwrap();
        let mut direct = Searcher::new(&ligand, None, &ff, search_box).unwrap();
        direct.pair_tables = None;
        let mut rng = StdRng::seed_from_u64(9);
        for _ in 0..20 {
            let conf = tabled.random_conf(&mut rng);
            let coords = tabled.coords_of(&conf);
            let a = tabled.intra(&coords, AUTHENTIC_CAP, None).unwrap();
            let b = direct.intra(&coords, AUTHENTIC_CAP, None).unwrap();
            assert!(
                (a - b).abs() < 1e-3 * (1.0 + b.abs()),
                "table {a} vs direct {b}"
            );
        }
    }

    /// The identity pose reproduces the input coordinates.
    #[test]
    fn test_identity_conf_is_input_pose() {
        let ligand = data("ligand.pdbqt");
        let ff = VinaForceField::default();
        let search_box = SearchBox::new(&Vector3::zeros(), &Vector3::repeat(20.0), 0.375);
        let searcher = Searcher::new(&ligand, None, &ff, search_box).unwrap();
        let coords = searcher.coords_of(&searcher.model.identity());
        for (atom, p) in ligand.atoms.iter().zip(&coords) {
            assert!((atom.coordinates - p).norm() < 1e-12);
        }
    }

    /// Out-of-box atoms are pushed back: the penalty gradient points outward.
    #[test]
    fn test_search_box_penalizes_outside() {
        let b = SearchBox::new(&Vector3::zeros(), &Vector3::repeat(10.0), 0.375);
        let (adjusted, miss, sign) = b.clamp(&Vector3::new(7.0, 0.0, -8.0));
        assert!(miss > 0.0);
        assert_eq!(sign, Vector3::new(1.0, 0.0, -1.0));
        assert!(b.contains(&adjusted, 1e-9));
    }

    #[test]
    fn test_container_keeps_distinct_low_poses() {
        let pose = |x: f64, e: f64| Pose {
            conf: Conf::identity(0, 0),
            e,
            heavy: vec![Vector3::new(x, 0.0, 0.0)],
        };
        let mut out = Vec::new();
        add_to_container(&mut out, pose(0.0, -1.0), 1.0, 2);
        add_to_container(&mut out, pose(0.5, -2.0), 1.0, 2); // same basin, better
        add_to_container(&mut out, pose(5.0, -0.5), 1.0, 2);
        add_to_container(&mut out, pose(9.0, -3.0), 1.0, 2); // full: replaces worst
        let energies: Vec<f64> = out.iter().map(|p| p.e).collect();
        assert_eq!(energies, vec![-3.0, -2.0]);
    }
}
