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
pub trait ForceField {
    /// Get the name of the forcefield
    fn name(&self) -> &'static str;

    /// Calculate the interaction energy between two atoms
    fn atom_pair_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError>;

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

    /// Calculate the interaction between an atom and a point charge
    fn atom_charge_energy(
        &self,
        atom: &Atom,
        charge: f64,
        position: &Vector3<f64>,
    ) -> Result<f64, ForceFieldError>;
}
