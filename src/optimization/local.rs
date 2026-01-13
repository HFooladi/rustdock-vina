//! Local optimization module using gradient-based methods
//!
//! This module implements BFGS-style local optimization for refining docking poses.

use nalgebra::{DVector, Unit, UnitQuaternion, Vector3};

use crate::forcefield::ForceField;
use crate::molecule::Molecule;
use crate::optimization::OptimizationError;

/// Parameters for local optimization
#[derive(Debug, Clone)]
pub struct LocalOptimizerParams {
    /// Step size for finite difference gradient computation
    pub gradient_step: f64,

    /// Initial step size for line search
    pub initial_step: f64,

    /// Convergence tolerance for energy
    pub energy_tolerance: f64,

    /// Convergence tolerance for gradient
    pub gradient_tolerance: f64,

    /// Maximum number of iterations
    pub max_iterations: usize,

    /// Armijo parameter for line search (c1)
    pub armijo_c1: f64,

    /// Wolfe parameter for line search (c2)
    pub wolfe_c2: f64,
}

impl Default for LocalOptimizerParams {
    fn default() -> Self {
        Self {
            gradient_step: 0.01,      // Angstroms for translation, radians for rotation
            initial_step: 0.1,        // Initial step size
            energy_tolerance: 1e-6,   // kcal/mol
            gradient_tolerance: 1e-4, // Gradient norm threshold
            max_iterations: 100,      // Max BFGS iterations
            armijo_c1: 1e-4,          // Armijo condition parameter
            wolfe_c2: 0.9,            // Wolfe condition parameter
        }
    }
}

/// Local optimizer using L-BFGS method
pub struct LocalOptimizer {
    pub params: LocalOptimizerParams,
}

impl Default for LocalOptimizer {
    fn default() -> Self {
        Self {
            params: LocalOptimizerParams::default(),
        }
    }
}

impl LocalOptimizer {
    /// Create a new local optimizer with default parameters
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a new local optimizer with custom parameters
    pub fn with_params(params: LocalOptimizerParams) -> Self {
        Self { params }
    }

    /// Perform local optimization on a molecule
    ///
    /// Optimizes translation, rotation, and torsion angles to minimize energy.
    pub fn minimize(
        &self,
        ligand: &mut Molecule,
        receptor: &Molecule,
        forcefield: &dyn ForceField,
    ) -> Result<f64, OptimizationError> {
        // Get number of degrees of freedom: 3 translation + 3 rotation + n_torsions
        let n_torsions = ligand.torsions.len();
        let _n_dof = 6 + n_torsions; // 3 trans + 3 rot + torsions

        // Get initial state
        let mut x = self.get_state_vector(ligand);
        let mut energy = self.calculate_energy(ligand, receptor, forcefield)?;

        // L-BFGS history
        let m = 10; // Number of corrections to store
        let mut s_history: Vec<DVector<f64>> = Vec::with_capacity(m);
        let mut y_history: Vec<DVector<f64>> = Vec::with_capacity(m);
        let mut rho_history: Vec<f64> = Vec::with_capacity(m);

        let mut prev_grad = self.compute_gradient(ligand, receptor, forcefield, &x)?;

        for iter in 0..self.params.max_iterations {
            // Compute gradient
            let grad = if iter == 0 {
                prev_grad.clone()
            } else {
                self.compute_gradient(ligand, receptor, forcefield, &x)?
            };

            // Check convergence
            let grad_norm = grad.norm();
            if grad_norm < self.params.gradient_tolerance {
                break;
            }

            // Compute search direction using L-BFGS two-loop recursion
            let direction = self.lbfgs_direction(&grad, &s_history, &y_history, &rho_history);

            // Line search
            let (step_size, new_energy) =
                self.line_search(ligand, receptor, forcefield, &x, &direction, energy, &grad)?;

            if step_size < 1e-10 {
                // Line search failed, stop optimization
                break;
            }

            // Update state
            let new_x = &x + step_size * &direction;

            // Update L-BFGS history
            let s = &new_x - &x;
            let new_grad = self.compute_gradient(ligand, receptor, forcefield, &new_x)?;
            let y = &new_grad - &grad;

            let sy = s.dot(&y);
            if sy > 1e-10 {
                // Only add to history if curvature condition is satisfied
                if s_history.len() >= m {
                    s_history.remove(0);
                    y_history.remove(0);
                    rho_history.remove(0);
                }
                s_history.push(s);
                y_history.push(y.clone());
                rho_history.push(1.0 / sy);
            }

            // Check energy convergence
            let energy_change = (new_energy - energy).abs();
            if energy_change < self.params.energy_tolerance {
                self.apply_state_vector(ligand, &new_x);
                energy = new_energy;
                break;
            }

            x = new_x;
            energy = new_energy;
            prev_grad = new_grad;

            // Apply state to molecule
            self.apply_state_vector(ligand, &x);
        }

        Ok(energy)
    }

    /// Get the current state as a vector [tx, ty, tz, rx, ry, rz, torsion1, torsion2, ...]
    fn get_state_vector(&self, molecule: &Molecule) -> DVector<f64> {
        let n_torsions = molecule.torsions.len();
        let mut state = DVector::zeros(6 + n_torsions);

        // Translation is relative (starts at 0)
        state[0] = 0.0;
        state[1] = 0.0;
        state[2] = 0.0;

        // Rotation is relative (starts at 0)
        state[3] = 0.0;
        state[4] = 0.0;
        state[5] = 0.0;

        // Torsion angles
        for (i, torsion) in molecule.torsions.iter().enumerate() {
            state[6 + i] = torsion.angle;
        }

        state
    }

    /// Apply state vector to molecule
    fn apply_state_vector(&self, molecule: &mut Molecule, state: &DVector<f64>) {
        // Apply translation
        let translation = Vector3::new(state[0], state[1], state[2]);
        for atom in &mut molecule.atoms {
            atom.coordinates += translation;
        }

        // Apply rotation (axis-angle representation)
        let rot_axis = Vector3::new(state[3], state[4], state[5]);
        let rot_angle = rot_axis.norm();

        if rot_angle > 1e-10 {
            let unit_axis = Unit::new_normalize(rot_axis);
            let rotation = UnitQuaternion::from_axis_angle(&unit_axis, rot_angle);

            // Rotate around center of mass
            if let Ok(center) = molecule.center() {
                for atom in &mut molecule.atoms {
                    let pos = atom.coordinates - center;
                    let rotated = rotation.transform_vector(&pos);
                    atom.coordinates = rotated + center;
                }
            }
        }

        // Apply torsion angles
        for (i, torsion) in molecule.torsions.iter_mut().enumerate() {
            torsion.angle = state[6 + i];
        }
    }

    /// Compute gradient using finite differences
    fn compute_gradient(
        &self,
        molecule: &Molecule,
        receptor: &Molecule,
        forcefield: &dyn ForceField,
        _state: &DVector<f64>,
    ) -> Result<DVector<f64>, OptimizationError> {
        let n_torsions = molecule.torsions.len();
        let n_dof = 6 + n_torsions;
        let mut grad = DVector::zeros(n_dof);

        let h = self.params.gradient_step;
        let base_energy = self.calculate_energy(molecule, receptor, forcefield)?;

        // Gradient for translation (x, y, z)
        for i in 0..3 {
            let mut mol_plus = molecule.clone();
            let mut delta = Vector3::zeros();
            delta[i] = h;
            for atom in &mut mol_plus.atoms {
                atom.coordinates += delta;
            }
            let energy_plus = self.calculate_energy(&mol_plus, receptor, forcefield)?;
            grad[i] = (energy_plus - base_energy) / h;
        }

        // Gradient for rotation (rx, ry, rz)
        for i in 0..3 {
            let mut mol_plus = molecule.clone();
            let mut axis = Vector3::zeros();
            axis[i] = 1.0;
            let unit_axis = Unit::new_normalize(axis);
            let rotation = UnitQuaternion::from_axis_angle(&unit_axis, h);

            if let Ok(center) = mol_plus.center() {
                for atom in &mut mol_plus.atoms {
                    let pos = atom.coordinates - center;
                    let rotated = rotation.transform_vector(&pos);
                    atom.coordinates = rotated + center;
                }
            }
            let energy_plus = self.calculate_energy(&mol_plus, receptor, forcefield)?;
            grad[3 + i] = (energy_plus - base_energy) / h;
        }

        // Gradient for torsions
        for i in 0..n_torsions {
            let mut mol_plus = molecule.clone();
            // Apply small torsion change
            if i < mol_plus.torsions.len() {
                mol_plus.torsions[i].angle += h;
                // Note: In a full implementation, we would also rotate the atoms
                // For now, this is a simplified version
            }
            let energy_plus = self.calculate_energy(&mol_plus, receptor, forcefield)?;
            grad[6 + i] = (energy_plus - base_energy) / h;
        }

        Ok(grad)
    }

    /// L-BFGS two-loop recursion to compute search direction
    fn lbfgs_direction(
        &self,
        grad: &DVector<f64>,
        s_history: &[DVector<f64>],
        y_history: &[DVector<f64>],
        rho_history: &[f64],
    ) -> DVector<f64> {
        if s_history.is_empty() {
            // If no history, use steepest descent
            return -grad.clone();
        }

        let k = s_history.len();
        let mut q = grad.clone();
        let mut alpha = vec![0.0; k];

        // First loop (backward)
        for i in (0..k).rev() {
            alpha[i] = rho_history[i] * s_history[i].dot(&q);
            q = &q - alpha[i] * &y_history[i];
        }

        // Initial Hessian approximation (scaled identity)
        let gamma =
            s_history[k - 1].dot(&y_history[k - 1]) / y_history[k - 1].dot(&y_history[k - 1]);
        let mut r = gamma * q;

        // Second loop (forward)
        for i in 0..k {
            let beta = rho_history[i] * y_history[i].dot(&r);
            r = &r + (alpha[i] - beta) * &s_history[i];
        }

        -r
    }

    /// Backtracking line search with Armijo condition
    fn line_search(
        &self,
        molecule: &Molecule,
        receptor: &Molecule,
        forcefield: &dyn ForceField,
        x: &DVector<f64>,
        direction: &DVector<f64>,
        current_energy: f64,
        grad: &DVector<f64>,
    ) -> Result<(f64, f64), OptimizationError> {
        let mut step = self.params.initial_step;
        let c1 = self.params.armijo_c1;
        let rho = 0.5; // Step reduction factor

        let directional_derivative = grad.dot(direction);

        // If direction is not a descent direction, return failure
        if directional_derivative >= 0.0 {
            return Ok((0.0, current_energy));
        }

        for _ in 0..20 {
            // Max line search iterations
            let new_x = x + step * direction;

            // Create molecule copy and apply new state
            let mut new_mol = molecule.clone();
            self.apply_state_vector(&mut new_mol, &new_x);

            let new_energy = self.calculate_energy(&new_mol, receptor, forcefield)?;

            // Armijo condition
            if new_energy <= current_energy + c1 * step * directional_derivative {
                return Ok((step, new_energy));
            }

            step *= rho;
        }

        // Line search failed
        Ok((0.0, current_energy))
    }

    /// Calculate energy for a molecule configuration
    fn calculate_energy(
        &self,
        ligand: &Molecule,
        receptor: &Molecule,
        forcefield: &dyn ForceField,
    ) -> Result<f64, OptimizationError> {
        let mut total_energy = 0.0;
        let cutoff = 8.0;

        // Intermolecular energy
        for lig_atom in &ligand.atoms {
            for rec_atom in &receptor.atoms {
                let distance = lig_atom.distance(rec_atom);
                if distance < cutoff && distance > 0.01 {
                    let pair_energy = forcefield.atom_pair_energy(lig_atom, rec_atom, distance)?;
                    total_energy += pair_energy;
                }
            }
        }

        // Note: Internal energy not included since it cancels with unbound energy in Vina formula

        Ok(total_energy)
    }
}
