//! Grid representation for efficient force field evaluation

use crate::atom::{Atom, AtomType};
use crate::forcefield::{ForceField, ForceFieldError};
use crate::molecule::Molecule;
use nalgebra::Vector3;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::path::Path;
use thiserror::Error;

/// Errors that can occur when working with grids
#[derive(Error, Debug)]
pub enum GridError {
    #[error("Invalid grid dimension: {0}")]
    InvalidDimension(String),

    #[error("Point {0:?} is outside grid bounds")]
    OutOfBounds(Vector3<f64>),

    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("Parse error: {0}")]
    ParseError(String),
}

/// A 3D grid for storing precomputed potentials
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Grid {
    /// Origin of the grid (minimum corner) in Angstroms
    pub origin: Vector3<f64>,

    /// Grid point spacing in Angstroms
    pub spacing: f64,

    /// Number of grid points in each dimension
    pub dimensions: Vector3<usize>,

    /// Atom type this grid represents
    pub atom_type: AtomType,

    /// Grid data - stored as a flattened 3D array
    pub data: Vec<f64>,
}

impl Grid {
    /// Create a new grid with the specified dimensions and initialize all values to zero
    pub fn new(
        origin: Vector3<f64>,
        spacing: f64,
        dimensions: Vector3<usize>,
        atom_type: AtomType,
    ) -> Result<Self, GridError> {
        if dimensions.x == 0 || dimensions.y == 0 || dimensions.z == 0 {
            return Err(GridError::InvalidDimension(format!("{:?}", dimensions)));
        }

        let total_points = dimensions.x * dimensions.y * dimensions.z;

        Ok(Self {
            origin,
            spacing,
            dimensions,
            atom_type,
            data: vec![0.0; total_points],
        })
    }

    /// Get the linear index for a 3D grid position
    fn get_index(&self, x: usize, y: usize, z: usize) -> Result<usize, GridError> {
        if x >= self.dimensions.x || y >= self.dimensions.y || z >= self.dimensions.z {
            return Err(GridError::OutOfBounds(Vector3::new(
                x as f64, y as f64, z as f64,
            )));
        }

        // Row-major order: (x * ny * nz) + (y * nz) + z
        Ok((x * self.dimensions.y * self.dimensions.z) + (y * self.dimensions.z) + z)
    }

    /// Get a value from the grid at integer coordinates
    pub fn get_value_at_indices(&self, x: usize, y: usize, z: usize) -> Result<f64, GridError> {
        let idx = self.get_index(x, y, z)?;
        Ok(self.data[idx])
    }

    /// Set a value in the grid at integer coordinates
    pub fn set_value_at_indices(
        &mut self,
        x: usize,
        y: usize,
        z: usize,
        value: f64,
    ) -> Result<(), GridError> {
        let idx = self.get_index(x, y, z)?;
        self.data[idx] = value;
        Ok(())
    }

    /// Convert real-space coordinates to grid indices
    pub fn real_to_grid(
        &self,
        position: &Vector3<f64>,
    ) -> Result<(usize, usize, usize), GridError> {
        // Calculate grid indices
        let x = ((position.x - self.origin.x) / self.spacing) as usize;
        let y = ((position.y - self.origin.y) / self.spacing) as usize;
        let z = ((position.z - self.origin.z) / self.spacing) as usize;

        // Check if within bounds
        if x >= self.dimensions.x || y >= self.dimensions.y || z >= self.dimensions.z {
            return Err(GridError::OutOfBounds(*position));
        }

        Ok((x, y, z))
    }

    /// Get a value from the grid at real-space coordinates with trilinear interpolation
    pub fn get_value(&self, position: &Vector3<f64>) -> Result<f64, GridError> {
        Ok(self.value_and_gradient(position)?.0)
    }

    /// Trilinearly interpolated value at `position` together with its
    /// gradient with respect to `position`.
    ///
    /// The gradient is the exact derivative of the interpolant, so it is
    /// piecewise constant along each axis within a cell and consistent with
    /// [`get_value`](Self::get_value).
    pub fn value_and_gradient(
        &self,
        position: &Vector3<f64>,
    ) -> Result<(f64, Vector3<f64>), GridError> {
        // Calculate grid indices as floating point
        let x_f = (position.x - self.origin.x) / self.spacing;
        let y_f = (position.y - self.origin.y) / self.spacing;
        let z_f = (position.z - self.origin.z) / self.spacing;

        // Check if within bounds (with a small margin for floating point errors)
        if x_f < -1e-6
            || y_f < -1e-6
            || z_f < -1e-6
            || x_f >= (self.dimensions.x as f64 - 1e-6)
            || y_f >= (self.dimensions.y as f64 - 1e-6)
            || z_f >= (self.dimensions.z as f64 - 1e-6)
        {
            return Err(GridError::OutOfBounds(*position));
        }

        // Get integer indices
        let x0 = x_f.max(0.0).floor() as usize;
        let y0 = y_f.max(0.0).floor() as usize;
        let z0 = z_f.max(0.0).floor() as usize;

        // Handle edge case at the maximum boundary
        let x1 = (x0 + 1).min(self.dimensions.x - 1);
        let y1 = (y0 + 1).min(self.dimensions.y - 1);
        let z1 = (z0 + 1).min(self.dimensions.z - 1);

        // Calculate interpolation factors
        let dx = x_f - x0 as f64;
        let dy = y_f - y0 as f64;
        let dz = z_f - z0 as f64;

        // Get the eight corner values
        let v000 = self.get_value_at_indices(x0, y0, z0)?;
        let v001 = self.get_value_at_indices(x0, y0, z1)?;
        let v010 = self.get_value_at_indices(x0, y1, z0)?;
        let v011 = self.get_value_at_indices(x0, y1, z1)?;
        let v100 = self.get_value_at_indices(x1, y0, z0)?;
        let v101 = self.get_value_at_indices(x1, y0, z1)?;
        let v110 = self.get_value_at_indices(x1, y1, z0)?;
        let v111 = self.get_value_at_indices(x1, y1, z1)?;

        // Perform trilinear interpolation
        let c00 = v000 * (1.0 - dx) + v100 * dx;
        let c01 = v001 * (1.0 - dx) + v101 * dx;
        let c10 = v010 * (1.0 - dx) + v110 * dx;
        let c11 = v011 * (1.0 - dx) + v111 * dx;

        let c0 = c00 * (1.0 - dy) + c10 * dy;
        let c1 = c01 * (1.0 - dy) + c11 * dy;

        let value = c0 * (1.0 - dz) + c1 * dz;

        // Partial derivatives of the interpolant in grid units, then scaled
        // to real space.
        let d_dz = c1 - c0;
        let d_dy = (c10 - c00) * (1.0 - dz) + (c11 - c01) * dz;
        let e00 = v100 - v000;
        let e01 = v101 - v001;
        let e10 = v110 - v010;
        let e11 = v111 - v011;
        let d_dx = (e00 * (1.0 - dy) + e10 * dy) * (1.0 - dz) + (e01 * (1.0 - dy) + e11 * dy) * dz;

        Ok((value, Vector3::new(d_dx, d_dy, d_dz) / self.spacing))
    }

    /// Vina's grid evaluation (`grid::evaluate`): the interpolant at
    /// `position` clamped into the grid, softened by Vina's `curl` with cap `cap`,
    /// plus `slope` times the distance outside the grid along each axis.
    ///
    /// Returns the value and its gradient. Along an axis where the point is
    /// outside, the gradient is the penalty's alone, pushing it back in.
    pub fn evaluate_clamped(
        &self,
        position: &Vector3<f64>,
        slope: f64,
        cap: f64,
    ) -> (f64, Vector3<f64>) {
        let dims = [self.dimensions.x, self.dimensions.y, self.dimensions.z];
        let mut s = [0.0; 3];
        let mut cell = [0usize; 3];
        let mut region = [0.0; 3];
        let mut miss = 0.0;
        for i in 0..3 {
            let si = (position[i] - self.origin[i]) / self.spacing;
            let last = (dims[i] - 1) as f64;
            if si < 0.0 {
                miss += -si;
                region[i] = -1.0;
                cell[i] = 0;
                s[i] = 0.0;
            } else if si >= last {
                miss += si - last;
                region[i] = 1.0;
                cell[i] = dims[i] - 2;
                s[i] = 1.0;
            } else {
                cell[i] = si as usize;
                s[i] = si - cell[i] as f64;
            }
        }
        let penalty = slope * miss * self.spacing;

        let (ny, nz) = (dims[1], dims[2]);
        let at = |x: usize, y: usize, z: usize| self.data[x * ny * nz + y * nz + z];
        let [x0, y0, z0] = cell;
        let (f000, f100) = (at(x0, y0, z0), at(x0 + 1, y0, z0));
        let (f010, f110) = (at(x0, y0 + 1, z0), at(x0 + 1, y0 + 1, z0));
        let (f001, f101) = (at(x0, y0, z0 + 1), at(x0 + 1, y0, z0 + 1));
        let (f011, f111) = (at(x0, y0 + 1, z0 + 1), at(x0 + 1, y0 + 1, z0 + 1));
        let [x, y, z] = s;
        let (mx, my, mz) = (1.0 - x, 1.0 - y, 1.0 - z);

        let mut f = f000 * mx * my * mz
            + f100 * x * my * mz
            + f010 * mx * y * mz
            + f110 * x * y * mz
            + f001 * mx * my * z
            + f101 * x * my * z
            + f011 * mx * y * z
            + f111 * x * y * z;
        let mut g = Vector3::new(
            (f100 - f000) * my * mz
                + (f110 - f010) * y * mz
                + (f101 - f001) * my * z
                + (f111 - f011) * y * z,
            (f010 - f000) * mx * mz
                + (f110 - f100) * x * mz
                + (f011 - f001) * mx * z
                + (f111 - f101) * x * z,
            (f001 - f000) * mx * my
                + (f101 - f100) * x * my
                + (f011 - f010) * mx * y
                + (f111 - f110) * x * y,
        );
        curl(&mut f, &mut g, cap);
        for i in 0..3 {
            let inside = if region[i] == 0.0 { g[i] } else { 0.0 };
            g[i] = inside / self.spacing + slope * region[i];
        }
        (f + penalty, g)
    }

    /// Load a grid from an AutoDock map file
    pub fn from_file<P: AsRef<Path>>(_path: P, _atom_type: AtomType) -> Result<Self, GridError> {
        // This would normally parse AutoDock map files
        // For now, return an error as this is not implemented
        Err(GridError::ParseError("Not implemented".to_string()))
    }

    /// Save grid to an AutoDock map file
    pub fn to_file<P: AsRef<Path>>(&self, _path: P) -> Result<(), GridError> {
        // This would normally save to AutoDock map format
        // For now, return an error as this is not implemented
        Err(GridError::ParseError("Not implemented".to_string()))
    }
}

/// What a ligand atom looks like to a type-gridable force field: its type
/// plus the bond-graph-derived `xs` flags. Two atoms with the same key score
/// identically against any receptor atom at any distance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct ProbeKey {
    atom_type: AtomType,
    is_hydrophobic: bool,
    is_hbond_donor: bool,
}

impl ProbeKey {
    fn of(atom: &Atom) -> Self {
        Self {
            atom_type: atom.atom_type,
            is_hydrophobic: atom.is_hydrophobic,
            is_hbond_donor: atom.is_hbond_donor,
        }
    }
}

/// Precomputed ligand-receptor interaction grids, one per ligand atom class.
///
/// This is how Vina makes the search affordable: each grid point holds the
/// summed pair energy of a probe atom against every receptor atom in range,
/// so scoring a pose costs one interpolation per ligand atom instead of a sum
/// over every receptor atom.
///
/// Atoms that wander outside the gridded region are scored exactly against
/// the stored receptor atoms, so leaving the grid costs speed, not accuracy.
#[derive(Debug, Clone)]
pub struct GridSet {
    grids: Vec<Grid>,
    /// Grid index for each ligand atom; `None` for atoms that never interact
    /// with the receptor (hydrogens and closure dummies).
    atom_grid: Vec<Option<usize>>,
    /// Receptor heavy atoms, kept for scoring atoms outside the grid.
    receptor: Vec<Atom>,
}

impl GridSet {
    /// Build grids covering the search box widened by `margin` on every side.
    ///
    /// Grids are built only for the atom classes present in `ligand`, and the
    /// set is tied to that ligand's atom order. Returns `Ok(None)` when the
    /// force field's pair energy cannot be gridded
    /// (see [`ForceField::grid_cutoff`]).
    pub fn build(
        ligand: &Molecule,
        receptor: &Molecule,
        forcefield: &dyn ForceField,
        center: &Vector3<f64>,
        box_size: &Vector3<f64>,
        margin: f64,
        spacing: f64,
    ) -> Result<Option<Self>, GridError> {
        let Some(cutoff) = forcefield.grid_cutoff() else {
            return Ok(None);
        };
        if !spacing.is_finite() || spacing <= 0.0 {
            return Err(GridError::InvalidDimension(format!(
                "grid spacing must be positive, got {spacing}"
            )));
        }

        // As Vina: a whole number of cells spanning the box, centred on it.
        let cells = (box_size + Vector3::repeat(2.0 * margin))
            .map(|extent| ((extent / spacing).ceil() as usize).max(1));
        let origin = center - cells.map(|n| n as f64 * spacing) * 0.5;
        let dimensions = cells.map(|n| n + 1);

        let mut probes: Vec<Atom> = Vec::new();
        let mut keys: Vec<ProbeKey> = Vec::new();
        let atom_grid = ligand
            .atoms
            .iter()
            .map(|atom| {
                // Vina's grids skip closure dummies as well as hydrogens.
                if atom.atom_type.is_hydrogen() || atom.atom_type.is_closure_dummy() {
                    return None;
                }
                let key = ProbeKey::of(atom);
                Some(keys.iter().position(|k| *k == key).unwrap_or_else(|| {
                    keys.push(key);
                    probes.push(atom.clone());
                    keys.len() - 1
                }))
            })
            .collect();

        let receptor: Vec<Atom> = receptor
            .atoms
            .iter()
            .filter(|a| !a.atom_type.is_hydrogen())
            .cloned()
            .collect();

        let grids = match probes.first() {
            Some(first) => {
                let template = Grid::new(origin, spacing, dimensions, first.atom_type)?;
                fill_grids(&template, &probes, &receptor, forcefield, cutoff)?
            }
            None => Vec::new(),
        };

        Ok(Some(Self {
            grids,
            atom_grid,
            receptor,
        }))
    }

    /// Number of grids, i.e. distinct interacting atom classes in the ligand.
    pub fn len(&self) -> usize {
        self.grids.len()
    }

    /// Whether the set holds no grids (a ligand with no heavy atoms).
    pub fn is_empty(&self) -> bool {
        self.grids.is_empty()
    }

    /// Vina's search energy (`cache::eval_deriv`) of ligand atoms at
    /// `coords`, adding each atom's gradient into `grads`.
    ///
    /// Every atom is evaluated on its grid with [`Grid::evaluate_clamped`]:
    /// an atom outside the grid is charged `slope` per Angstrom, which is how
    /// Vina confines the whole ligand to the search box. `cap` softens each
    /// atom's positive energy (see [`curl`]).
    pub(crate) fn search_energy(
        &self,
        coords: &[Vector3<f64>],
        slope: f64,
        cap: f64,
        grads: &mut [Vector3<f64>],
    ) -> f64 {
        let mut total = 0.0;
        for ((p, grid), grad) in coords.iter().zip(&self.atom_grid).zip(grads) {
            let Some(g) = grid else { continue };
            let (e, d) = self.grids[*g].evaluate_clamped(p, slope, cap);
            total += e;
            *grad += d;
        }
        total
    }

    /// Ligand-receptor interaction energy of `ligand`, which must have the
    /// same atoms in the same order as the ligand the set was built for.
    pub fn energy(
        &self,
        ligand: &Molecule,
        forcefield: &dyn ForceField,
    ) -> Result<f64, ForceFieldError> {
        let mut total = 0.0;
        for (atom, grid) in ligand.atoms.iter().zip(&self.atom_grid) {
            let Some(g) = grid else { continue };
            total += match self.grids[*g].get_value(&atom.coordinates) {
                Ok(value) => value,
                Err(_) => pairwise_energy(atom, &self.receptor, forcefield)?,
            };
        }
        Ok(total)
    }

    /// As [`energy`](Self::energy), additionally *adding* each ligand atom's
    /// energy gradient `dE/dr_i` into `grads` (one entry per ligand atom).
    pub fn energy_and_gradients(
        &self,
        ligand: &Molecule,
        forcefield: &dyn ForceField,
        grads: &mut [Vector3<f64>],
    ) -> Result<f64, ForceFieldError> {
        let mut total = 0.0;
        for ((atom, grid), grad) in ligand.atoms.iter().zip(&self.atom_grid).zip(grads) {
            let Some(g) = grid else { continue };
            let (energy, gradient) = match self.grids[*g].value_and_gradient(&atom.coordinates) {
                Ok(found) => found,
                Err(_) => pairwise_energy_and_gradient(atom, &self.receptor, forcefield)?,
            };
            total += energy;
            *grad += gradient;
        }
        Ok(total)
    }
}

/// Vina's `curl`: soften a positive energy `e` with cap `v` to
/// `e * v / (v + e)`, scaling its derivative by the square of that factor.
///
/// It keeps clashes from dominating the search: a pose with one bad contact
/// is penalized but not driven off a cliff. Negative energies pass unchanged.
pub(crate) fn curl(e: &mut f64, deriv: &mut Vector3<f64>, v: f64) {
    if *e > 0.0 && v.is_finite() {
        let factor = if v < f64::EPSILON { 0.0 } else { v / (v + *e) };
        *e *= factor;
        *deriv *= factor * factor;
    }
}

/// Exact energy of one atom against a list of receptor atoms, and its
/// gradient with respect to the atom's position.
pub(crate) fn pairwise_energy_and_gradient(
    atom: &Atom,
    receptor: &[Atom],
    forcefield: &dyn ForceField,
) -> Result<(f64, Vector3<f64>), ForceFieldError> {
    let mut energy = 0.0;
    let mut gradient = Vector3::zeros();
    for rec_atom in receptor {
        let delta = atom.coordinates - rec_atom.coordinates;
        let distance = delta.norm();
        if distance > 0.01 {
            energy += forcefield.atom_pair_energy(atom, rec_atom, distance)?;
            gradient +=
                delta * (forcefield.atom_pair_derivative(atom, rec_atom, distance)? / distance);
        }
    }
    Ok((energy, gradient))
}

/// Exact energy of one atom against a list of receptor atoms.
fn pairwise_energy(
    atom: &Atom,
    receptor: &[Atom],
    forcefield: &dyn ForceField,
) -> Result<f64, ForceFieldError> {
    let mut total = 0.0;
    for rec_atom in receptor {
        let distance = atom.distance(rec_atom);
        if distance > 0.01 {
            total += forcefield.atom_pair_energy(atom, rec_atom, distance)?;
        }
    }
    Ok(total)
}

/// Resolution of the radial pair-energy tables, in Angstroms.
///
/// Linear interpolation at this step is accurate to ~1e-6 kcal/mol on the
/// smooth terms and ~2e-4 at the kinks of piecewise-linear ones, far below the
/// trilinear error of the grid itself.
const TABLE_STEP: f64 = 1.0 / 1024.0;

/// Pair energy against distance for every (receptor class, ligand class)
/// pair, sampled every [`TABLE_STEP`] up to the cutoff.
///
/// Laid out as `[receptor class][sample][ligand class]`, so one receptor atom
/// at one distance reads its contribution to every ligand class from two
/// contiguous runs.
struct RadialTables {
    data: Vec<f64>,
    samples: usize,
    classes: usize,
}

impl RadialTables {
    fn build(
        probes: &[Atom],
        receptor_classes: &[Atom],
        forcefield: &dyn ForceField,
        cutoff: f64,
    ) -> Result<Self, ForceFieldError> {
        // One sample past the cutoff, so the upper interpolation neighbour of
        // any distance <= cutoff exists.
        let samples = (cutoff / TABLE_STEP).floor() as usize + 2;
        let classes = probes.len();
        let mut data = Vec::with_capacity(receptor_classes.len() * samples * classes);
        for rec in receptor_classes {
            for k in 0..samples {
                // Clamped to the range the grid fill evaluates, so samples
                // outside it carry the edge value rather than a cutoff zero.
                let r = (k as f64 * TABLE_STEP).clamp(MIN_DISTANCE, cutoff);
                for probe in probes {
                    data.push(forcefield.atom_pair_energy(probe, rec, r)?);
                }
            }
        }
        Ok(Self {
            data,
            samples,
            classes,
        })
    }

    /// Add the energy of receptor class `rec` at distance `r` to `out`, one
    /// entry per ligand class.
    #[inline]
    fn accumulate(&self, rec: usize, r: f64, out: &mut [f64]) {
        let x = r / TABLE_STEP;
        let k = (x as usize).min(self.samples - 2);
        let t = x - k as f64;
        let base = (rec * self.samples + k) * self.classes;
        let (lo, hi) = self.data[base..base + 2 * self.classes].split_at(self.classes);
        for ((v, a), b) in out.iter_mut().zip(lo).zip(hi) {
            *v += a + t * (b - a);
        }
    }
}

/// Pairs closer than this are skipped, as in exact scoring.
const MIN_DISTANCE: f64 = 0.01;

/// Fill one grid per probe with the probe's energy at every grid point, in a
/// single sweep over the receptor.
///
/// Parallel over x-slabs. Each row keeps only receptor atoms within `cutoff`
/// of the row's axis, sorted by z, and each point then reads the window of
/// that list within `cutoff` in z. Every pair costs one square root and a
/// table lookup per probe, with no force-field call.
fn fill_grids(
    template: &Grid,
    probes: &[Atom],
    receptor: &[Atom],
    forcefield: &dyn ForceField,
    cutoff: f64,
) -> Result<Vec<Grid>, GridError> {
    let table_error =
        |e: ForceFieldError| GridError::ParseError(format!("grid energy evaluation failed: {e}"));

    // Receptor atoms by class, dropping any too far from the grid to matter.
    let lo = template.origin - Vector3::repeat(cutoff);
    let hi = template.origin
        + template.dimensions.map(|n| (n - 1) as f64) * template.spacing
        + Vector3::repeat(cutoff);
    let mut class_keys: Vec<ProbeKey> = Vec::new();
    let mut class_atoms: Vec<Atom> = Vec::new();
    let mut sites: Vec<(Vector3<f64>, usize)> = Vec::new();
    for atom in receptor {
        let p = atom.coordinates;
        if (0..3).any(|i| p[i] < lo[i] || p[i] > hi[i]) {
            continue;
        }
        let key = ProbeKey::of(atom);
        let class = class_keys
            .iter()
            .position(|k| *k == key)
            .unwrap_or_else(|| {
                class_keys.push(key);
                class_atoms.push(atom.clone());
                class_keys.len() - 1
            });
        sites.push((p, class));
    }
    let tables =
        RadialTables::build(probes, &class_atoms, forcefield, cutoff).map_err(table_error)?;

    let n = probes.len();
    let (ny, nz) = (template.dimensions.y, template.dimensions.z);
    let (origin, spacing) = (template.origin, template.spacing);
    let cutoff2 = cutoff * cutoff;
    let min2 = MIN_DISTANCE * MIN_DISTANCE;

    // Interleaved `[x][y][z][probe]`, so a point's values are contiguous.
    let mut values = vec![0.0; template.data.len() * n];
    values
        .par_chunks_mut(ny * nz * n)
        .enumerate()
        .for_each(|(ix, slab)| {
            let x = origin.x + ix as f64 * spacing;
            let near_x: Vec<(Vector3<f64>, usize)> = sites
                .iter()
                .copied()
                .filter(|(p, _)| (p.x - x).abs() <= cutoff)
                .collect();
            // (z, squared distance from the row's axis, class)
            let mut row_sites: Vec<(f64, f64, usize)> = Vec::with_capacity(near_x.len());

            for (iy, row) in slab.chunks_mut(nz * n).enumerate() {
                let y = origin.y + iy as f64 * spacing;
                row_sites.clear();
                for (p, class) in &near_x {
                    let axial2 = (p.x - x).powi(2) + (p.y - y).powi(2);
                    if axial2 <= cutoff2 {
                        row_sites.push((p.z, axial2, *class));
                    }
                }
                row_sites.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));

                let (mut first, mut last) = (0, 0);
                for (iz, point) in row.chunks_mut(n).enumerate() {
                    let z = origin.z + iz as f64 * spacing;
                    while first < row_sites.len() && row_sites[first].0 < z - cutoff {
                        first += 1;
                    }
                    while last < row_sites.len() && row_sites[last].0 <= z + cutoff {
                        last += 1;
                    }
                    for &(sz, axial2, class) in &row_sites[first..last] {
                        let r2 = axial2 + (sz - z).powi(2);
                        if r2 > min2 && r2 <= cutoff2 {
                            tables.accumulate(class, r2.sqrt(), point);
                        }
                    }
                }
            }
        });

    probes
        .iter()
        .enumerate()
        .map(|(i, probe)| {
            let mut grid = template.clone();
            grid.atom_type = probe.atom_type;
            for (v, point) in grid.data.iter_mut().zip(values.chunks_exact(n)) {
                *v = point[i];
            }
            Ok(grid)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::{ad4::AD4ForceField, vina::VinaForceField};

    fn atom(atom_type: AtomType, x: f64, y: f64, z: f64) -> Atom {
        Atom::new(
            atom_type,
            Vector3::new(x, y, z),
            "TEST".to_string(),
            1,
            "TST".to_string(),
            1,
            'A',
            0.0,
        )
    }

    fn receptor() -> Molecule {
        let mut m = Molecule::new("rec");
        m.add_atom(atom(AtomType::Carbon, 0.0, 0.0, 0.0));
        m.add_atom(atom(AtomType::OxygenH, 3.0, 1.0, -1.0));
        m.add_atom(atom(AtomType::NitrogenH, -2.0, 2.5, 1.5));
        m
    }

    fn ligand_at(positions: &[Vector3<f64>]) -> Molecule {
        let mut m = Molecule::new("lig");
        for (i, p) in positions.iter().enumerate() {
            let t = if i % 2 == 0 {
                AtomType::Carbon
            } else {
                AtomType::OxygenH
            };
            m.add_atom(atom(t, p.x, p.y, p.z));
        }
        m
    }

    fn exact(ligand: &Molecule, receptor: &Molecule, ff: &dyn ForceField) -> f64 {
        ligand
            .atoms
            .iter()
            .map(|a| pairwise_energy(a, &receptor.atoms, ff).unwrap())
            .sum()
    }

    fn build(ligand: &Molecule, spacing: f64) -> GridSet {
        GridSet::build(
            ligand,
            &receptor(),
            &VinaForceField::default(),
            &Vector3::zeros(),
            &Vector3::repeat(12.0),
            0.0,
            spacing,
        )
        .unwrap()
        .expect("Vina is gridable")
    }

    #[test]
    fn test_grid_matches_pairwise_on_grid_points() {
        // Box [-6, 6] with spacing 0.5 puts points on every half-angstrom.
        let ligand = ligand_at(&[Vector3::new(4.0, 0.5, 0.0), Vector3::new(1.5, 3.0, -2.0)]);
        let grids = build(&ligand, 0.5);
        let ff = VinaForceField::default();

        // On grid points only the radial tables' interpolation separates the
        // two.
        let grid_energy = grids.energy(&ligand, &ff).unwrap();
        assert!((grid_energy - exact(&ligand, &receptor(), &ff)).abs() < 1e-3);
    }

    #[test]
    fn test_every_grid_point_matches_pairwise() {
        // A receptor dense and varied enough to exercise every term,
        // including the hydrophobic and H-bond ramps of flagged atoms.
        let types = [
            AtomType::Carbon,
            AtomType::OxygenH,
            AtomType::NitrogenH,
            AtomType::Nitrogen,
            AtomType::SulfurH,
        ];
        let mut rec = Molecule::new("rec");
        let mut state = 7u64;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        for i in 0..60 {
            let mut a = atom(
                types[i % types.len()],
                next() * 16.0 - 8.0,
                next() * 16.0 - 8.0,
                next() * 16.0 - 8.0,
            );
            a.is_hydrophobic = a.atom_type == AtomType::Carbon && i % 2 == 0;
            a.is_hbond_donor = a.atom_type == AtomType::NitrogenH && i % 3 == 0;
            rec.add_atom(a);
        }
        let mut ligand = ligand_at(&[
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
        ]);
        ligand.atoms[0].is_hydrophobic = true;
        ligand.atoms[2].is_hydrophobic = false;
        ligand.atoms[3].is_hbond_donor = true;
        let ff = VinaForceField::default();
        let grids = GridSet::build(
            &ligand,
            &rec,
            &ff,
            &Vector3::zeros(),
            &Vector3::repeat(10.0),
            0.0,
            0.5,
        )
        .unwrap()
        .unwrap();
        assert_eq!(grids.len(), 4);

        let mut worst: f64 = 0.0;
        for (atom, g) in ligand.atoms.iter().zip(&grids.atom_grid) {
            let grid = &grids.grids[g.unwrap()];
            let d = grid.dimensions;
            for ix in 0..d.x {
                for iy in 0..d.y {
                    for iz in 0..d.z {
                        let mut probe = atom.clone();
                        probe.coordinates = grid.origin
                            + Vector3::new(ix as f64, iy as f64, iz as f64) * grid.spacing;
                        let expected = pairwise_energy(&probe, &rec.atoms, &ff).unwrap();
                        let got = grid.get_value_at_indices(ix, iy, iz).unwrap();
                        worst = worst.max((got - expected).abs() / expected.abs().max(1.0));
                    }
                }
            }
        }
        assert!(worst < 1e-3, "worst relative grid error {worst}");
    }

    #[test]
    fn test_grid_interpolation_is_close_off_grid() {
        let ligand = ligand_at(&[Vector3::new(4.1, 0.7, 0.3), Vector3::new(1.37, 3.21, -2.11)]);
        let grids = build(&ligand, 0.375);
        let ff = VinaForceField::default();

        let grid_energy = grids.energy(&ligand, &ff).unwrap();
        let exact_energy = exact(&ligand, &receptor(), &ff);
        assert!(
            (grid_energy - exact_energy).abs() < 0.01,
            "grid {grid_energy} vs exact {exact_energy}"
        );
    }

    #[test]
    fn test_atoms_outside_grid_are_scored_exactly() {
        // 7.5 is outside the [-6, 6] box but within the 8 A cutoff of the
        // receptor oxygen, so a silent zero would be wrong.
        let outside = ligand_at(&[Vector3::new(7.5, 0.0, 0.0)]);
        let grids = build(&outside, 0.5);
        let ff = VinaForceField::default();

        let grid_energy = grids.energy(&outside, &ff).unwrap();
        let exact_energy = exact(&outside, &receptor(), &ff);
        assert!(exact_energy != 0.0, "test point must feel the receptor");
        assert!((grid_energy - exact_energy).abs() < 1e-12);
    }

    #[test]
    fn test_grid_gradient_matches_finite_difference() {
        let ligand = ligand_at(&[Vector3::new(4.1, 0.7, 0.3)]);
        let grids = build(&ligand, 0.375);
        let grid = &grids.grids[0];

        let h = 1e-6;
        for p in [
            Vector3::new(4.1, 0.7, 0.3),
            Vector3::new(1.37, 3.21, -2.11),
            Vector3::new(-2.9, -1.05, 4.4),
        ] {
            let (_, analytic) = grid.value_and_gradient(&p).unwrap();
            for axis in 0..3 {
                let mut step = Vector3::zeros();
                step[axis] = h;
                let numeric = (grid.get_value(&(p + step)).unwrap()
                    - grid.get_value(&(p - step)).unwrap())
                    / (2.0 * h);
                assert!(
                    (numeric - analytic[axis]).abs() < 1e-5,
                    "axis {axis} at {p:?}: analytic {} vs numeric {numeric}",
                    analytic[axis]
                );
            }
        }
    }

    #[test]
    fn test_out_of_grid_gradient_is_exact() {
        let ligand = ligand_at(&[Vector3::new(7.5, 0.0, 0.0)]);
        let grids = build(&ligand, 0.5);
        let ff = VinaForceField::default();
        let mut grads = vec![Vector3::zeros(); 1];
        grids
            .energy_and_gradients(&ligand, &ff, &mut grads)
            .unwrap();

        let (_, expected) =
            pairwise_energy_and_gradient(&ligand.atoms[0], &receptor().atoms, &ff).unwrap();
        assert!((grads[0] - expected).norm() < 1e-12);
        assert!(expected.norm() > 0.0, "test point must feel the receptor");
    }

    #[test]
    fn test_one_grid_per_atom_class() {
        let ligand = ligand_at(&[
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
        ]);
        assert_eq!(build(&ligand, 1.0).len(), 2);
    }

    #[test]
    fn test_charge_dependent_forcefield_is_not_gridded() {
        let ligand = ligand_at(&[Vector3::new(1.0, 0.0, 0.0)]);
        let grids = GridSet::build(
            &ligand,
            &receptor(),
            &AD4ForceField::default(),
            &Vector3::zeros(),
            &Vector3::repeat(12.0),
            0.0,
            0.5,
        )
        .unwrap();
        assert!(grids.is_none());
    }
}
