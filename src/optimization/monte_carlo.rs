//! Monte Carlo search algorithm for molecular docking

use nalgebra::{Vector3, UnitQuaternion, Unit};
use rand::prelude::*;
use std::f64::consts::PI;

use crate::forcefield::ForceField;
use crate::molecule::Molecule;
use crate::optimization::{DockingResult, OptimizationError, Optimizer};

/// Parameters for Monte Carlo search
#[derive(Debug, Clone)]
pub struct MonteCarloParams {
    /// Temperature parameter for Metropolis criterion
    pub temperature: f64,
    
    /// Probability of performing a translation move
    pub prob_translation: f64,
    
    /// Probability of performing a rotation move
    pub prob_rotation: f64,
    
    /// Probability of performing a torsion move
    pub prob_torsion: f64,
    
    /// Maximum translation step (in Angstroms)
    pub max_translation: f64,
    
    /// Maximum rotation step (in radians)
    pub max_rotation: f64,
    
    /// Maximum torsion step (in radians)
    pub max_torsion: f64,
    
    /// Number of steps without improvement before stopping
    pub steps_without_improvement: usize,
}

impl Default for MonteCarloParams {
    fn default() -> Self {
        Self {
            temperature: 0.3,
            prob_translation: 0.5,
            prob_rotation: 0.3,
            prob_torsion: 0.2,
            max_translation: 2.0,
            max_rotation: PI / 6.0, // 30 degrees
            max_torsion: PI / 6.0,  // 30 degrees
            steps_without_improvement: 1000,
        }
    }
}

/// Implementation of Monte Carlo search for molecular docking
#[derive(Debug, Clone)]
pub struct MonteCarlo {
    pub params: MonteCarloParams,
}

impl Default for MonteCarlo {
    fn default() -> Self {
        Self {
            params: MonteCarloParams::default(),
        }
    }
}

impl MonteCarlo {
    /// Create a new Monte Carlo optimizer with default parameters
    pub fn new() -> Self {
        Self::default()
    }
    
    /// Create a new Monte Carlo optimizer with custom parameters
    pub fn with_params(params: MonteCarloParams) -> Self {
        Self { params }
    }
    
    /// Generate a random translation within the search box
    fn random_translation(&self, rng: &mut ThreadRng) -> Vector3<f64> {
        let dx = (rng.gen::<f64>() - 0.5) * 2.0 * self.params.max_translation;
        let dy = (rng.gen::<f64>() - 0.5) * 2.0 * self.params.max_translation;
        let dz = (rng.gen::<f64>() - 0.5) * 2.0 * self.params.max_translation;
        
        Vector3::new(dx, dy, dz)
    }
    
    /// Generate a random rotation
    fn random_rotation(&self, rng: &mut ThreadRng) -> UnitQuaternion<f64> {
        let angle = (rng.gen::<f64>() - 0.5) * 2.0 * self.params.max_rotation;
        let axis = Vector3::new(
            rng.gen::<f64>() - 0.5,
            rng.gen::<f64>() - 0.5,
            rng.gen::<f64>() - 0.5,
        );
        
        if axis.norm() > 0.0 {
            let unit_axis = Unit::new_normalize(axis);
            UnitQuaternion::from_axis_angle(&unit_axis, angle)
        } else {
            // If axis is zero, return identity quaternion
            UnitQuaternion::identity()
        }
    }
    
    /// Apply a translation to the molecule
    fn translate_molecule(&self, molecule: &mut Molecule, translation: &Vector3<f64>) {
        for atom in &mut molecule.atoms {
            atom.coordinates += translation;
        }
    }
    
    /// Apply a rotation to the molecule around its center
    fn rotate_molecule(&self, molecule: &mut Molecule, rotation: &UnitQuaternion<f64>) {
        // Calculate center of molecule
        let center = match molecule.center() {
            Ok(c) => c,
            Err(_) => return, // No atoms in molecule
        };
        
        // Apply rotation around center
        for atom in &mut molecule.atoms {
            // Move to origin
            let pos = atom.coordinates - center;
            
            // Apply rotation
            let rotated = rotation.transform_vector(&pos);
            
            // Move back
            atom.coordinates = rotated + center;
        }
    }
    
    /// Modify a single torsion angle
    fn modify_torsion(&self, molecule: &mut Molecule, torsion_idx: usize, delta_angle: f64) -> Result<(), OptimizationError> {
        if torsion_idx >= molecule.torsions.len() {
            return Err(OptimizationError::InvalidInitialState);
        }
        
        let torsion = &mut molecule.torsions[torsion_idx];
        let bond_idx = torsion.bond_idx;
        
        if bond_idx >= molecule.bonds.len() {
            return Err(OptimizationError::InvalidInitialState);
        }
        
        let bond = &molecule.bonds[bond_idx];
        let atom1_idx = bond.atom1_idx;
        let atom2_idx = bond.atom2_idx;
        
        if atom1_idx >= molecule.atoms.len() || atom2_idx >= molecule.atoms.len() {
            return Err(OptimizationError::InvalidInitialState);
        }
        
        // Get positions of the atoms that define the bond
        let atom1_pos = molecule.atoms[atom1_idx].coordinates;
        let atom2_pos = molecule.atoms[atom2_idx].coordinates;
        
        // Calculate the axis of rotation
        let axis = Unit::new_normalize(atom2_pos - atom1_pos);
        
        // Create rotation quaternion
        let rotation = UnitQuaternion::from_axis_angle(&axis, delta_angle);
        
        // Apply rotation to all moving atoms
        for &atom_idx in &torsion.moving_atoms {
            if atom_idx >= molecule.atoms.len() {
                return Err(OptimizationError::InvalidInitialState);
            }
            
            // Move to origin relative to atom1
            let pos = molecule.atoms[atom_idx].coordinates - atom1_pos;
            
            // Apply rotation
            let rotated = rotation.transform_vector(&pos);
            
            // Move back
            molecule.atoms[atom_idx].coordinates = rotated + atom1_pos;
        }
        
        // Update torsion angle
        torsion.angle += delta_angle;
        
        // Normalize angle to [-PI, PI]
        while torsion.angle > PI {
            torsion.angle -= 2.0 * PI;
        }
        while torsion.angle < -PI {
            torsion.angle += 2.0 * PI;
        }
        
        Ok(())
    }
    
    /// Check if a molecule is within the search box
    fn is_in_box(&self, molecule: &Molecule, center: &Vector3<f64>, box_size: &Vector3<f64>) -> bool {
        let half_size = box_size * 0.5;
        let min_corner = center - half_size;
        let max_corner = center + half_size;
        
        for atom in &molecule.atoms {
            let pos = &atom.coordinates;
            
            if pos.x < min_corner.x || pos.x > max_corner.x ||
               pos.y < min_corner.y || pos.y > max_corner.y ||
               pos.z < min_corner.z || pos.z > max_corner.z {
                return false;
            }
        }
        
        true
    }
    
    /// Calculate energy using the provided forcefield
    fn calculate_energy(&self, molecule: &Molecule, forcefield: &dyn ForceField) -> Result<f64, OptimizationError> {
        // This is a placeholder - in a real implementation, this would calculate
        // the interaction energy between the molecule and the receptor
        
        // For now, just return the internal energy of the molecule
        let energy = molecule.calculate_internal_energy(8.0);
        
        Ok(energy)
    }
}

impl Optimizer for MonteCarlo {
    fn optimize(
        &self,
        molecule: &Molecule,
        forcefield: &dyn ForceField,
        center: Vector3<f64>,
        box_size: Vector3<f64>,
        max_steps: usize,
    ) -> Result<DockingResult, OptimizationError> {
        let mut rng = thread_rng();
        let mut current_molecule = molecule.clone();
        
        // Calculate initial energy
        let mut current_energy = self.calculate_energy(&current_molecule, forcefield)?;
        let mut best_energy = current_energy;
        let mut best_molecule = current_molecule.clone();
        
        let mut steps_since_improvement = 0;
        
        for step in 0..max_steps {
            // Create a copy of the current state
            let mut new_molecule = current_molecule.clone();
            
            // Choose a type of move
            let move_type = rng.gen::<f64>();
            
            if move_type < self.params.prob_translation {
                // Translation
                let translation = self.random_translation(&mut rng);
                self.translate_molecule(&mut new_molecule, &translation);
            } else if move_type < self.params.prob_translation + self.params.prob_rotation {
                // Rotation
                let rotation = self.random_rotation(&mut rng);
                self.rotate_molecule(&mut new_molecule, &rotation);
            } else {
                // Torsion
                if !new_molecule.torsions.is_empty() {
                    // Choose a random torsion
                    let torsion_idx = rng.gen_range(0..new_molecule.torsions.len());
                    let delta_angle = (rng.gen::<f64>() - 0.5) * 2.0 * self.params.max_torsion;
                    
                    // Apply torsion
                    if let Err(e) = self.modify_torsion(&mut new_molecule, torsion_idx, delta_angle) {
                        // Skip this step if torsion modification fails
                        continue;
                    }
                }
            }
            
            // Check if the molecule is still within the search box
            if !self.is_in_box(&new_molecule, &center, &box_size) {
                continue;
            }
            
            // Calculate new energy
            let new_energy = self.calculate_energy(&new_molecule, forcefield)?;
            
            // Decide whether to accept the move (Metropolis criterion)
            let accept = if new_energy <= current_energy {
                true
            } else {
                let boltzmann_factor = ((current_energy - new_energy) / self.params.temperature).exp();
                rng.gen::<f64>() < boltzmann_factor
            };
            
            if accept {
                current_molecule = new_molecule;
                current_energy = new_energy;
                
                // Check if this is a new best solution
                if new_energy < best_energy {
                    best_energy = new_energy;
                    best_molecule = current_molecule.clone();
                    steps_since_improvement = 0;
                } else {
                    steps_since_improvement += 1;
                }
            } else {
                steps_since_improvement += 1;
            }
            
            // Check termination criteria
            if steps_since_improvement >= self.params.steps_without_improvement {
                break;
            }
        }
        
        // Create energy components for the result
        let energy_components = vec![
            ("Total".to_string(), best_energy),
            // In a real implementation, this would include individual components
            // like vdW, electrostatic, etc.
        ];
        
        Ok(DockingResult {
            energy: best_energy,
            molecule: best_molecule,
            energy_components,
            rmsd: None,
        })
    }
    
    fn generate_poses(
        &self,
        molecule: &Molecule,
        forcefield: &dyn ForceField,
        center: Vector3<f64>,
        box_size: Vector3<f64>,
        num_poses: usize,
        max_steps: usize,
    ) -> Result<Vec<DockingResult>, OptimizationError> {
        let mut results = Vec::with_capacity(num_poses);
        
        // Generate poses
        for _ in 0..num_poses {
            let result = self.optimize(molecule, forcefield, center, box_size, max_steps)?;
            results.push(result);
        }
        
        // Sort by energy
        results.sort_by(|a, b| a.energy.partial_cmp(&b.energy).unwrap());
        
        // Calculate RMSD to best pose
        if !results.is_empty() {
            let (first, rest) = results.split_at_mut(1);
            let best_molecule = &first[0].molecule;
            
            for result in rest {
                let rmsd = calculate_rmsd(&result.molecule, best_molecule);
                result.rmsd = Some(rmsd);
            }
        }
        
        Ok(results)
    }
}

/// Calculate root mean square deviation between two molecules
fn calculate_rmsd(mol1: &Molecule, mol2: &Molecule) -> f64 {
    if mol1.atoms.len() != mol2.atoms.len() {
        return f64::MAX; // Different number of atoms
    }
    
    let mut sum_sq_diff = 0.0;
    
    for i in 0..mol1.atoms.len() {
        let pos1 = &mol1.atoms[i].coordinates;
        let pos2 = &mol2.atoms[i].coordinates;
        
        let diff = pos1 - pos2;
        sum_sq_diff += diff.norm_squared();
    }
    
    (sum_sq_diff / mol1.atoms.len() as f64).sqrt()
}