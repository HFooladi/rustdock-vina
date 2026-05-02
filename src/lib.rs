//! RustDock-Vina: A Rust implementation of the AutoDock Vina molecular docking algorithm.
//!
//! This library provides molecular docking simulations using the AutoDock Vina
//! scoring function and Monte Carlo optimization with local BFGS refinement.
//!
//! # Quick start
//!
//! Parse two PDBQT files and compute the intermolecular interaction energy
//! using the Vina force field:
//!
//! ```no_run
//! use rustdock_vina::{parse_pdbqt, ForceField, Result, VinaForceField};
//!
//! fn score_pose() -> Result<()> {
//!     let receptor = parse_pdbqt("receptor.pdbqt")?;
//!     let ligand = parse_pdbqt("ligand.pdbqt")?;
//!     let forcefield = VinaForceField::default();
//!
//!     let cutoff = 8.0;
//!     let mut energy = 0.0;
//!     for lig_atom in &ligand.atoms {
//!         for rec_atom in &receptor.atoms {
//!             let d = lig_atom.distance(rec_atom);
//!             if d < cutoff && d > 0.01 {
//!                 energy += forcefield.atom_pair_energy(lig_atom, rec_atom, d)?;
//!             }
//!         }
//!     }
//!     println!("Intermolecular energy: {:.3} kcal/mol", energy);
//!     Ok(())
//! }
//! ```
//!
//! # Error handling
//!
//! Each submodule defines its own concrete error type (e.g. [`IoError`],
//! [`ForceFieldError`]), and the crate-level [`enum@Error`] wraps all of them
//! via `#[from]`. This lets callers either let `?` aggregate failures into
//! [`Result<T>`], or pattern-match on the specific cause:
//!
//! ```
//! use rustdock_vina::{Error, Result, parse_pdbqt};
//!
//! fn load() -> Result<()> {
//!     let _ = parse_pdbqt("missing.pdbqt")?; // IoError -> Error via From
//!     Ok(())
//! }
//!
//! match load() {
//!     Ok(()) => {}
//!     Err(Error::Io(e)) => eprintln!("I/O failure: {e}"),
//!     Err(other) => eprintln!("Other failure: {other}"),
//! }
//! ```

pub mod atom;
pub mod forcefield;
pub mod grid;
pub mod io;
pub mod math;
pub mod molecule;
pub mod optimization;
pub mod scoring;
pub mod utils;

use thiserror::Error;

pub use atom::{Atom, AtomType};
pub use forcefield::ad4::{AD4ForceField, AD4Params};
pub use forcefield::vina::{VinaForceField, VinaParams};
pub use forcefield::{ForceField, ForceFieldError};
pub use grid::{Grid, GridError};
pub use io::{parse_pdbqt, write_docking_results, write_pdbqt, IoError};
pub use molecule::{Bond, Molecule, MoleculeError, Torsion};
pub use optimization::monte_carlo::{MonteCarlo, MonteCarloParams};
pub use optimization::{DockingResult, OptimizationError, Optimizer};
pub use scoring::vina_score::VinaScore;

/// Version of the library, sourced from `Cargo.toml`.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// All errors that can be returned from this crate's public API.
///
/// Each variant wraps a module-specific error from a submodule, allowing
/// callers to either treat failures uniformly via [`Display`](std::fmt::Display)
/// or pattern-match on the specific cause.
#[derive(Debug, Error)]
pub enum Error {
    #[error(transparent)]
    Io(#[from] IoError),

    #[error(transparent)]
    Molecule(#[from] MoleculeError),

    #[error(transparent)]
    ForceField(#[from] ForceFieldError),

    #[error(transparent)]
    Optimization(#[from] OptimizationError),

    #[error(transparent)]
    Grid(#[from] GridError),
}

/// Convenience alias for `Result<T, Error>` used throughout the public API.
pub type Result<T> = std::result::Result<T, Error>;
