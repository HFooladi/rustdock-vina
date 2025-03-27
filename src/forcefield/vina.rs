//! Implementation of the Vina scoring function

use crate::atom::{Atom, AtomType};
use crate::forcefield::{ForceField, ForceFieldError};
use nalgebra::Vector3;

/// Parameters for Vina forcefield
#[derive(Debug, Clone)]
pub struct VinaParams {
    // Weights for each component of the scoring function
    pub weight_gauss1: f64,
    pub weight_gauss2: f64,
    pub weight_repulsion: f64,
    pub weight_hydrophobic: f64,
    pub weight_hydrogen: f64,
    pub weight_rot: f64,
    
    // Gaussian function parameters
    pub gaussian_offset: f64,
    pub gaussian_width: f64,
    
    // Hydrogen bond parameters
    pub hydrogen_bond_dist_cutoff: f64,
    pub hydrogen_bond_angle_cutoff: f64,
}

impl Default for VinaParams {
    fn default() -> Self {
        // Default parameters from the Vina paper
        Self {
            weight_gauss1: -0.0356,
            weight_gauss2: -0.00516,
            weight_repulsion: 0.840,
            weight_hydrophobic: -0.0351,
            weight_hydrogen: -0.587,
            weight_rot: 0.05846,
            
            gaussian_offset: 0.0,
            gaussian_width: 0.5,
            
            hydrogen_bond_dist_cutoff: 4.0,
            hydrogen_bond_angle_cutoff: 120.0,
        }
    }
}

/// Implementation of the Vina scoring function
#[derive(Debug, Clone)]
pub struct VinaForceField {
    pub params: VinaParams,
}

impl Default for VinaForceField {
    fn default() -> Self {
        Self {
            params: VinaParams::default(),
        }
    }
}

impl ForceField for VinaForceField {
    fn name(&self) -> &'static str {
        "Vina"
    }
    
    fn atom_pair_energy(&self, atom1: &Atom, atom2: &Atom, distance: f64) -> Result<f64, ForceFieldError> {
        // Early return if distance is too large
        if distance > 8.0 {
            return Ok(0.0);
        }
        
        // Calculate all interaction components
        let vdw = self.vdw_energy(atom1, atom2, distance)?;
        let hbond = self.hbond_energy(atom1, atom2, distance)?;
        let hydrophobic = self.hydrophobic_energy(atom1, atom2, distance)?;
        let elec = self.electrostatic_energy(atom1, atom2, distance)?;
        let desol = self.desolvation_energy(atom1, atom2, distance)?;
        
        // Sum all components
        Ok(vdw + hbond + hydrophobic + elec + desol)
    }
    
    fn vdw_energy(&self, atom1: &Atom, atom2: &Atom, distance: f64) -> Result<f64, ForceFieldError> {
        // Skip if atoms are too close
        if distance < 0.01 {
            return Ok(0.0);
        }
        
        // Calculate optimal distance (sum of vdW radii)
        let r1 = atom1.atom_type.radius();
        let r2 = atom2.atom_type.radius();
        let optimal_dist = r1 + r2;
        
        // Gaussian attractive term
        let gauss1 = self.params.weight_gauss1 * 
            (-(distance - optimal_dist - self.params.gaussian_offset).powi(2) / 
             (2.0 * self.params.gaussian_width.powi(2))).exp();
        
        // Repulsive term (only when atoms are too close)
        let repulsion = if distance < optimal_dist {
            self.params.weight_repulsion * (optimal_dist - distance).powi(2)
        } else {
            0.0
        };
        
        Ok(gauss1 + repulsion)
    }
    
    fn electrostatic_energy(&self, atom1: &Atom, atom2: &Atom, distance: f64) -> Result<f64, ForceFieldError> {
        // Simplified electrostatic interaction with distance-dependent dielectric
        // In the original Vina, this is handled differently
        let e = 332.0 * atom1.charge * atom2.charge / (distance * distance);
        Ok(e)
    }
    
    fn hbond_energy(&self, atom1: &Atom, atom2: &Atom, distance: f64) -> Result<f64, ForceFieldError> {
        // Check if we have a hydrogen bond acceptor and donor
        let (donor, acceptor) = if atom1.is_h_bond_donor() && atom2.is_h_bond_acceptor() {
            (atom1, atom2)
        } else if atom2.is_h_bond_donor() && atom1.is_h_bond_acceptor() {
            (atom2, atom1)
        } else {
            // No hydrogen bond possible
            return Ok(0.0);
        };
        
        // Check distance criteria
        if distance > self.params.hydrogen_bond_dist_cutoff {
            return Ok(0.0);
        }
        
        // In a real implementation, we would check the angle here as well
        // For now, we'll use a simple distance-based model
        
        // Calculate hydrogen bond strength based on distance
        let optimal_dist = 1.9; // Typical H-bond distance
        let strength = if distance <= optimal_dist {
            1.0
        } else {
            let normalized_dist = (distance - optimal_dist) / 
                                 (self.params.hydrogen_bond_dist_cutoff - optimal_dist);
            1.0 - normalized_dist
        };
        
        Ok(self.params.weight_hydrogen * strength)
    }
    
    fn hydrophobic_energy(&self, atom1: &Atom, atom2: &Atom, distance: f64) -> Result<f64, ForceFieldError> {
        // Check if both atoms are hydrophobic
        let is_hydrophobic1 = matches!(atom1.atom_type, AtomType::Carbon);
        let is_hydrophobic2 = matches!(atom2.atom_type, AtomType::Carbon);
        
        if !is_hydrophobic1 || !is_hydrophobic2 {
            return Ok(0.0);
        }
        
        // Distance cutoffs for hydrophobic interactions
        let r_max = 4.5;
        let r_min = 0.5;
        
        if distance > r_max || distance < r_min {
            return Ok(0.0);
        }
        
        // Linear interpolation between cutoffs
        let hydrophobic_factor = if distance <= r_min {
            1.0
        } else if distance >= r_max {
            0.0
        } else {
            (r_max - distance) / (r_max - r_min)
        };
        
        Ok(self.params.weight_hydrophobic * hydrophobic_factor)
    }
    
    fn desolvation_energy(&self, atom1: &Atom, atom2: &Atom, distance: f64) -> Result<f64, ForceFieldError> {
        // Simplified desolvation model
        // In a real implementation, this would be more complex
        
        // Default solvation parameters
        let sol1 = match atom1.atom_type {
            AtomType::Carbon => 0.6,
            AtomType::Oxygen | AtomType::OxygenH => -1.0,
            AtomType::Nitrogen | AtomType::NitrogenH => -1.0,
            AtomType::Sulfur | AtomType::SulfurH => 0.6,
            _ => 0.0,
        };
        
        let sol2 = match atom2.atom_type {
            AtomType::Carbon => 0.6,
            AtomType::Oxygen | AtomType::OxygenH => -1.0,
            AtomType::Nitrogen | AtomType::NitrogenH => -1.0,
            AtomType::Sulfur | AtomType::SulfurH => 0.6,
            _ => 0.0,
        };
        
        // Calculate desolvation energy
        let desolv_cutoff = 8.0;
        
        if distance > desolv_cutoff {
            return Ok(0.0);
        }
        
        // Desolvation factor
        let factor = 1.0 - (distance / desolv_cutoff).powi(2);
        let energy = factor * sol1 * sol2;
        
        Ok(energy)
    }
    
    fn atom_charge_energy(&self, atom: &Atom, charge: f64, position: &Vector3<f64>) -> Result<f64, ForceFieldError> {
        let distance = (atom.coordinates - position).norm();
        
        // Early return if distance is too small or too large
        if distance < 0.01 || distance > 8.0 {
            return Ok(0.0);
        }
        
        // Coulomb's law with distance-dependent dielectric
        let energy = 332.0 * atom.charge * charge / (distance * distance);
        
        Ok(energy)
    }
}

// Implement additional methods specific to the Vina forcefield
impl VinaForceField {
    /// Calculate the conformational entropy penalty for rotatable bonds
    pub fn rotatable_bond_penalty(&self, num_rotatable_bonds: usize) -> f64 {
        self.params.weight_rot * num_rotatable_bonds as f64
    }
}