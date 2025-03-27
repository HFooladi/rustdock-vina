//! Molecule representation and related functionality

use crate::atom::{Atom, AtomType};
use nalgebra::Vector3;
use std::collections::HashMap;
use thiserror::Error;

/// Errors that can occur when working with molecules
#[derive(Error, Debug)]
pub enum MoleculeError {
    #[error("Invalid bond: atoms {0} and {1} not found")]
    InvalidBond(u32, u32),
    
    #[error("Invalid atom index: {0}")]
    InvalidAtomIndex(usize),
    
    #[error("No atoms in molecule")]
    EmptyMolecule,
}

/// Represents a chemical bond between two atoms
#[derive(Debug, Clone)]
pub struct Bond {
    /// Index of the first atom
    pub atom1_idx: usize,
    
    /// Index of the second atom
    pub atom2_idx: usize,
    
    /// Is this bond rotatable?
    pub rotatable: bool,
}

/// Represents a molecular torsion (rotatable bond and the atoms affected by rotation)
#[derive(Debug, Clone)]
pub struct Torsion {
    /// Index of the rotatable bond
    pub bond_idx: usize,
    
    /// Indices of atoms that move when this bond rotates
    pub moving_atoms: Vec<usize>,
    
    /// Current torsion angle (in radians)
    pub angle: f64,
}

/// Represents a molecule (ligand or receptor)
#[derive(Debug, Clone)]
pub struct Molecule {
    /// Name of the molecule
    pub name: String,
    
    /// List of atoms in the molecule
    pub atoms: Vec<Atom>,
    
    /// List of bonds between atoms
    pub bonds: Vec<Bond>,
    
    /// List of rotatable bonds (torsions)
    pub torsions: Vec<Torsion>,
    
    /// Is this molecule rigid or has flexible parts?
    pub is_rigid: bool,
}

impl Molecule {
    /// Create a new empty molecule
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            atoms: Vec::new(),
            bonds: Vec::new(),
            torsions: Vec::new(),
            is_rigid: true,
        }
    }
    
    /// Add an atom to the molecule
    pub fn add_atom(&mut self, atom: Atom) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(atom);
        idx
    }
    
    /// Add a bond between two atoms
    pub fn add_bond(&mut self, atom1_idx: usize, atom2_idx: usize, rotatable: bool) -> Result<usize, MoleculeError> {
        if atom1_idx >= self.atoms.len() || atom2_idx >= self.atoms.len() {
            return Err(MoleculeError::InvalidBond(atom1_idx as u32, atom2_idx as u32));
        }
        
        let bond = Bond {
            atom1_idx,
            atom2_idx,
            rotatable,
        };
        
        let idx = self.bonds.len();
        self.bonds.push(bond);
        
        if rotatable {
            self.is_rigid = false;
        }
        
        Ok(idx)
    }
    
    /// Get the center of the molecule
    pub fn center(&self) -> Result<Vector3<f64>, MoleculeError> {
        if self.atoms.is_empty() {
            return Err(MoleculeError::EmptyMolecule);
        }
        
        let sum = self.atoms.iter().fold(Vector3::zeros(), |acc, atom| {
            acc + atom.coordinates
        });
        
        Ok(sum / self.atoms.len() as f64)
    }
    
    /// Get the bounding box of the molecule
    pub fn bounding_box(&self) -> Result<(Vector3<f64>, Vector3<f64>), MoleculeError> {
        if self.atoms.is_empty() {
            return Err(MoleculeError::EmptyMolecule);
        }
        
        let mut min = Vector3::new(f64::MAX, f64::MAX, f64::MAX);
        let mut max = Vector3::new(f64::MIN, f64::MIN, f64::MIN);
        
        for atom in &self.atoms {
            // Update minimum corner
            min.x = min.x.min(atom.coordinates.x);
            min.y = min.y.min(atom.coordinates.y);
            min.z = min.z.min(atom.coordinates.z);
            
            // Update maximum corner
            max.x = max.x.max(atom.coordinates.x);
            max.y = max.y.max(atom.coordinates.y);
            max.z = max.z.max(atom.coordinates.z);
        }
        
        Ok((min, max))
    }
    
    /// Calculate atomic pairwise interactions within the molecule
    pub fn calculate_internal_energy(&self, max_distance: f64) -> f64 {
        let mut energy = 0.0;
        
        for i in 0..self.atoms.len() {
            for j in (i+1)..self.atoms.len() {
                let atom1 = &self.atoms[i];
                let atom2 = &self.atoms[j];
                
                // Skip atoms that are within 3 bonds of each other (1-4 interactions and closer)
                // This would normally be determined from the bond graph
                
                // Calculate distance
                let distance = atom1.distance(atom2);
                
                // Skip if beyond cutoff
                if distance > max_distance {
                    continue;
                }
                
                // In a real implementation, we would calculate van der Waals, 
                // electrostatic, and other interactions here
                // This is just a placeholder
                energy += 1.0 / distance;
            }
        }
        
        energy
    }
    
    /// Count atoms by type
    pub fn count_atom_types(&self) -> HashMap<AtomType, usize> {
        let mut counts = HashMap::new();
        
        for atom in &self.atoms {
            *counts.entry(atom.atom_type).or_insert(0) += 1;
        }
        
        counts
    }
}