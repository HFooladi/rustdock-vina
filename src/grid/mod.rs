//! Grid representation for efficient force field evaluation

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use std::path::Path;
use thiserror::Error;
use crate::atom::AtomType;

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
                x as f64,
                y as f64,
                z as f64,
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
    pub fn set_value_at_indices(&mut self, x: usize, y: usize, z: usize, value: f64) -> Result<(), GridError> {
        let idx = self.get_index(x, y, z)?;
        self.data[idx] = value;
        Ok(())
    }
    
    /// Convert real-space coordinates to grid indices
    pub fn real_to_grid(&self, position: &Vector3<f64>) -> Result<(usize, usize, usize), GridError> {
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
        if x_f < -1e-6 || y_f < -1e-6 || z_f < -1e-6 ||
           x_f >= (self.dimensions.x as f64 - 1e-6) || 
           y_f >= (self.dimensions.y as f64 - 1e-6) || 
           z_f >= (self.dimensions.z as f64 - 1e-6) {
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