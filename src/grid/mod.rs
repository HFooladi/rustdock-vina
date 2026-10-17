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
        let x0 = x_f.floor() as usize;
        let y0 = y_f.floor() as usize;
        let z0 = z_f.floor() as usize;

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

        let result = c0 * (1.0 - dz) + c1 * dz;

        Ok(result)
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
    /// (hydrogens).
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

        let half = box_size * 0.5 + Vector3::repeat(margin);
        let origin = center - half;
        let dimensions = (half * 2.0).map(|extent| (extent / spacing).ceil() as usize + 1);

        let mut probes: Vec<Atom> = Vec::new();
        let mut keys: Vec<ProbeKey> = Vec::new();
        let atom_grid = ligand
            .atoms
            .iter()
            .map(|atom| {
                if atom.atom_type.is_hydrogen() {
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

        let grids = probes
            .iter()
            .map(|probe| {
                let mut grid = Grid::new(origin, spacing, dimensions, probe.atom_type)?;
                fill_grid(&mut grid, probe, &receptor, forcefield, cutoff)?;
                Ok(grid)
            })
            .collect::<Result<Vec<_>, GridError>>()?;

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

/// Fill every point of `grid` with the energy of `probe` placed there.
///
/// Parallel over x-slabs. Each slab, then each row, narrows the receptor list
/// to atoms within `cutoff` along that axis, so the innermost loop only sees
/// a thin tube of candidates rather than the whole receptor.
fn fill_grid(
    grid: &mut Grid,
    probe: &Atom,
    receptor: &[Atom],
    forcefield: &dyn ForceField,
    cutoff: f64,
) -> Result<(), GridError> {
    let (ny, nz) = (grid.dimensions.y, grid.dimensions.z);
    let (origin, spacing) = (grid.origin, grid.spacing);

    grid.data
        .par_chunks_mut(ny * nz)
        .enumerate()
        .try_for_each(|(ix, slab)| {
            let x = origin.x + ix as f64 * spacing;
            let near_x: Vec<&Atom> = receptor
                .iter()
                .filter(|a| (a.coordinates.x - x).abs() <= cutoff)
                .collect();
            let mut probe = probe.clone();

            for (iy, row) in slab.chunks_mut(nz).enumerate() {
                let y = origin.y + iy as f64 * spacing;
                let near_xy: Vec<&Atom> = near_x
                    .iter()
                    .copied()
                    .filter(|a| (a.coordinates.y - y).abs() <= cutoff)
                    .collect();

                for (iz, value) in row.iter_mut().enumerate() {
                    probe.coordinates = Vector3::new(x, y, origin.z + iz as f64 * spacing);
                    let mut total = 0.0;
                    for rec_atom in &near_xy {
                        let distance = probe.distance(rec_atom);
                        if distance > 0.01 && distance <= cutoff {
                            total += forcefield.atom_pair_energy(&probe, rec_atom, distance)?;
                        }
                    }
                    *value = total;
                }
            }
            Ok::<(), ForceFieldError>(())
        })
        .map_err(|e| GridError::ParseError(format!("grid energy evaluation failed: {e}")))
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

        let grid_energy = grids.energy(&ligand, &ff).unwrap();
        assert!((grid_energy - exact(&ligand, &receptor(), &ff)).abs() < 1e-9);
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
        let ligand = ligand_at(&[Vector3::new(7.5, 0.0, 0.0), Vector3::new(1.5, 3.0, -2.0)]);
        let grids = build(&ligand, 0.5);
        let ff = VinaForceField::default();

        let grid_energy = grids.energy(&ligand, &ff).unwrap();
        assert!((grid_energy - exact(&ligand, &receptor(), &ff)).abs() < 1e-9);
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
