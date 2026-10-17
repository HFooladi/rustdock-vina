//! Forcefield implementations for molecular interactions

pub mod ad4;
pub mod vina;
use crate::atom::{Atom, AtomType};
use nalgebra::Vector3;
/// pub mod vinardo;
use thiserror::Error;

/// Errors that can occur in forcefields
#[derive(Error, Debug)]
pub enum ForceFieldError {
    #[error("Invalid atoms for interaction")]
    InvalidAtoms,

    #[error("Unsupported atom type: {0:?}")]
    UnsupportedAtomType(AtomType),

    #[error("Unknown error: {0}")]
    Unknown(String),
}

/// Trait representing a forcefield that can calculate interaction energies
/// The Send + Sync bounds enable parallel docking with rayon
pub trait ForceField: Send + Sync {
    /// Get the name of the forcefield
    fn name(&self) -> &'static str;

    /// Calculate the interaction energy between two atoms
    fn atom_pair_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError>;

    /// Derivative of [`atom_pair_energy`](Self::atom_pair_energy) with
    /// respect to the distance, `dE/dr`.
    ///
    /// The default is a central difference, which is correct for any force
    /// field; implementations with a closed form should override it.
    fn atom_pair_derivative(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError> {
        const H: f64 = 1e-5;
        let plus = self.atom_pair_energy(atom1, atom2, distance + H)?;
        let minus = self.atom_pair_energy(atom1, atom2, (distance - H).max(0.0))?;
        Ok((plus - minus) / (distance + H - (distance - H).max(0.0)))
    }

    /// Calculate van der Waals interaction energy
    fn vdw_energy(&self, atom1: &Atom, atom2: &Atom, distance: f64)
        -> Result<f64, ForceFieldError>;

    /// Calculate electrostatic interaction energy
    fn electrostatic_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError>;

    /// Calculate hydrogen bond energy
    fn hbond_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError>;

    /// Calculate desolvation energy
    fn desolvation_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError>;

    /// Calculate hydrophobic interaction energy
    fn hydrophobic_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError>;

    /// Calculate the interaction between an atom and a point charge
    fn atom_charge_energy(
        &self,
        atom: &Atom,
        charge: f64,
        position: &Vector3<f64>,
    ) -> Result<f64, ForceFieldError>;

    /// Apply torsional normalization to convert intermolecular energy to binding affinity
    /// This is used by Vina scoring: affinity = energy / (1 + weight_rot * N_rot)
    /// Default implementation returns the energy unchanged (for forcefields without this feature)
    fn apply_torsional_normalization(
        &self,
        intermolecular_energy: f64,
        _num_rotatable_bonds: usize,
    ) -> f64 {
        intermolecular_energy
    }

    /// Whether ligand-receptor energy can be precomputed on per-type grids,
    /// and over what range.
    ///
    /// Return `Some(cutoff)` only if [`atom_pair_energy`](Self::atom_pair_energy)
    /// depends on nothing but the two atoms' types, their `is_hydrophobic` and
    /// `is_hbond_donor` flags and the distance, is zero for hydrogens, and is
    /// zero beyond `cutoff`. Charge-dependent force fields must return `None`
    /// (the default), which keeps the search on exact pairwise sums.
    fn grid_cutoff(&self) -> Option<f64> {
        None
    }
}
