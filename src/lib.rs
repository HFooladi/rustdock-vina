//! RustDock-Vina: A Rust implementation of the AutoDock Vina molecular docking algorithm
//!
//! This library provides functionality for molecular docking simulations using the
//! AutoDock Vina scoring function and optimization algorithms.

pub mod atom;
pub mod forcefield;
pub mod grid;
pub mod io;
pub mod math;
pub mod molecule;
pub mod optimization;
pub mod scoring;
pub mod utils;

// Re-export commonly used types and functions
pub use atom::Atom;
pub use molecule::Molecule;
pub use scoring::vina_score::VinaScore;

/// Version of the library
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Result type for operations that can fail
pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
