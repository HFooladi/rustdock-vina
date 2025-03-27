//! Optimization algorithms for molecular docking

pub mod monte_carlo;

use crate::forcefield::ForceField;
use crate::molecule::Molecule;
use nalgebra::Vector3;
use thiserror::Error;

/// Errors that can occur during optimization
#[derive(Error, Debug)]
pub enum OptimizationError {
    #[error("Invalid initial state")]
    InvalidInitialState,

    #[error("Optimization failed to converge")]
    FailedToConverge,

    #[error("Maximum iterations exceeded")]
    MaxIterationsExceeded,

    #[error("ForceField error: {0}")]
    ForceFieldError(#[from] crate::forcefield::ForceFieldError),

    #[error("Unknown error: {0}")]
    Unknown(String),
}

/// Represents a docking result
#[derive(Debug, Clone)]
pub struct DockingResult {
    /// Final energy
    pub energy: f64,

    /// Docked ligand conformation
    pub molecule: Molecule,

    /// Components of the energy (vdW, electrostatic, etc.)
    pub energy_components: Vec<(String, f64)>,

    /// RMSD from the best mode (if multiple modes are generated)
    pub rmsd: Option<f64>,
}

/// Trait for optimization algorithms
pub trait Optimizer {
    /// Optimize the position, orientation and conformation of a molecule
    fn optimize(
        &self,
        molecule: &Molecule,
        forcefield: &dyn ForceField,
        center: Vector3<f64>,
        box_size: Vector3<f64>,
        max_steps: usize,
    ) -> Result<DockingResult, OptimizationError>;

    /// Generate multiple docking poses
    fn generate_poses(
        &self,
        molecule: &Molecule,
        forcefield: &dyn ForceField,
        center: Vector3<f64>,
        box_size: Vector3<f64>,
        num_poses: usize,
        max_steps: usize,
    ) -> Result<Vec<DockingResult>, OptimizationError>;
}
