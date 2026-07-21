//! Local optimization module using gradient-based methods
//!
//! This module implements BFGS-style local optimization for refining docking poses.

use nalgebra::{DVector, Unit, UnitQuaternion, Vector3};

use crate::forcefield::ForceField;
use crate::molecule::Molecule;
use crate::optimization::monte_carlo::{intermolecular_energy, intramolecular_energy};
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
#[derive(Default)]
pub struct LocalOptimizer {
    pub params: LocalOptimizerParams,
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

    /// Perform local optimization on a molecule.
    ///
    /// Optimizes translation, rotation, and torsion angles to minimize energy.
    /// Equivalent to [`minimize_in_box`](Self::minimize_in_box) with no box
    /// constraint.
    pub fn minimize(
        &self,
        ligand: &mut Molecule,
        receptor: &Molecule,
        forcefield: &dyn ForceField,
    ) -> Result<f64, OptimizationError> {
        self.minimize_in_box(ligand, receptor, forcefield, None)
    }

    /// Local optimization confined to a search box.
    ///
    /// `box_bounds` is `(center, size)`. Any trial pose that puts an atom
    /// outside the box is rejected. Without this the minimizer can walk the
    /// ligand clean out of the search volume and report the near-zero energy of
    /// a pose that is no longer touching the receptor.
    ///
    /// The state vector is `[dx, dy, dz, rx, ry, rz, dtor1, ...]` and every
    /// component is a *displacement from the pose the ligand had on entry*.
    /// Trial poses are always built by applying the whole state to an untouched
    /// copy of that entry pose, which keeps the mapping from state to
    /// coordinates a pure function — a prerequisite for the finite-difference
    /// gradients and the L-BFGS curvature updates to mean anything.
    pub fn minimize_in_box(
        &self,
        ligand: &mut Molecule,
        receptor: &Molecule,
        forcefield: &dyn ForceField,
        box_bounds: Option<(&Vector3<f64>, &Vector3<f64>)>,
    ) -> Result<f64, OptimizationError> {
        let pairs = ligand.intramolecular_pairs();
        self.minimize_in_box_with_intra(ligand, receptor, forcefield, box_bounds, &pairs)
    }

    /// As [`minimize_in_box`](Self::minimize_in_box), but reusing a
    /// precomputed list of intramolecular pairs.
    ///
    /// The pair list depends only on topology, so a caller running many
    /// minimizations on the same ligand should compute it once.
    pub fn minimize_in_box_with_intra(
        &self,
        ligand: &mut Molecule,
        receptor: &Molecule,
        forcefield: &dyn ForceField,
        box_bounds: Option<(&Vector3<f64>, &Vector3<f64>)>,
        intra_pairs: &[(usize, usize)],
    ) -> Result<f64, OptimizationError> {
        let base = ligand.clone();
        let n_dof = 6 + base.torsions.len();

        let mut x = DVector::zeros(n_dof);
        let mut energy = self.calculate_energy(&base, receptor, forcefield, intra_pairs)?;

        // L-BFGS history
        let m = 10; // Number of corrections to store
        let mut s_history: Vec<DVector<f64>> = Vec::with_capacity(m);
        let mut y_history: Vec<DVector<f64>> = Vec::with_capacity(m);
        let mut rho_history: Vec<f64> = Vec::with_capacity(m);

        let mut grad = self.compute_gradient(&base, receptor, forcefield, &x, intra_pairs)?;

        for _ in 0..self.params.max_iterations {
            if grad.norm() < self.params.gradient_tolerance {
                break;
            }

            let direction = self.lbfgs_direction(&grad, &s_history, &y_history, &rho_history);

            let (step_size, new_energy) = self.line_search(
                &base,
                receptor,
                forcefield,
                &x,
                &direction,
                energy,
                &grad,
                box_bounds,
                intra_pairs,
            )?;

            if step_size < 1e-10 {
                // Line search failed, stop optimization
                break;
            }

            let new_x = &x + step_size * &direction;
            let new_grad =
                self.compute_gradient(&base, receptor, forcefield, &new_x, intra_pairs)?;

            // L-BFGS curvature pair. Both gradients are taken at genuinely
            // different states, so `y` is non-zero and the history actually
            // fills; evaluating both around the same coordinates would silently
            // reduce this to steepest descent.
            let s = &new_x - &x;
            let y = &new_grad - &grad;
            let sy = s.dot(&y);
            if sy > 1e-10 {
                if s_history.len() >= m {
                    s_history.remove(0);
                    y_history.remove(0);
                    rho_history.remove(0);
                }
                s_history.push(s);
                y_history.push(y);
                rho_history.push(1.0 / sy);
            }

            let converged = (new_energy - energy).abs() < self.params.energy_tolerance;
            x = new_x;
            energy = new_energy;
            grad = new_grad;

            if converged {
                break;
            }
        }

        *ligand = self.pose_at(&base, &x);
        Ok(energy)
    }

    /// Build the conformation reached by applying state `x` to `base`.
    ///
    /// Torsions are applied first so that the rigid-body rotation and
    /// translation act on the final internal geometry.
    fn pose_at(&self, base: &Molecule, x: &DVector<f64>) -> Molecule {
        let mut mol = base.clone();

        for i in 0..base.torsions.len() {
            let delta = x[6 + i];
            if delta != 0.0 {
                // A torsion that cannot be rotated (malformed tree) is skipped
                // rather than aborting the whole minimization.
                let _ = mol.rotate_torsion(i, delta);
            }
        }

        let rot_axis = Vector3::new(x[3], x[4], x[5]);
        let rot_angle = rot_axis.norm();
        if rot_angle > 1e-10 {
            let rotation =
                UnitQuaternion::from_axis_angle(&Unit::new_normalize(rot_axis), rot_angle);
            mol.rotate_about_center(&rotation);
        }

        mol.translate(&Vector3::new(x[0], x[1], x[2]));
        mol
    }

    /// Compute gradient using forward finite differences in state space.
    fn compute_gradient(
        &self,
        base: &Molecule,
        receptor: &Molecule,
        forcefield: &dyn ForceField,
        state: &DVector<f64>,
        intra_pairs: &[(usize, usize)],
    ) -> Result<DVector<f64>, OptimizationError> {
        let n_dof = state.len();
        let mut grad = DVector::zeros(n_dof);

        let h = self.params.gradient_step;
        let base_energy = self.calculate_energy(
            &self.pose_at(base, state),
            receptor,
            forcefield,
            intra_pairs,
        )?;

        for i in 0..n_dof {
            let mut perturbed = state.clone();
            perturbed[i] += h;
            let energy_plus = self.calculate_energy(
                &self.pose_at(base, &perturbed),
                receptor,
                forcefield,
                intra_pairs,
            )?;
            grad[i] = (energy_plus - base_energy) / h;
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
    #[allow(clippy::too_many_arguments)]
    fn line_search(
        &self,
        base: &Molecule,
        receptor: &Molecule,
        forcefield: &dyn ForceField,
        x: &DVector<f64>,
        direction: &DVector<f64>,
        current_energy: f64,
        grad: &DVector<f64>,
        box_bounds: Option<(&Vector3<f64>, &Vector3<f64>)>,
        intra_pairs: &[(usize, usize)],
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
            let candidate = self.pose_at(base, &new_x);

            // Shrink the step rather than accepting a pose that leaves the box.
            let in_box = box_bounds
                .map(|(center, size)| candidate.center_within_box(center, size))
                .unwrap_or(true);

            if in_box {
                let new_energy =
                    self.calculate_energy(&candidate, receptor, forcefield, intra_pairs)?;

                // Armijo condition
                if new_energy <= current_energy + c1 * step * directional_derivative {
                    return Ok((step, new_energy));
                }
            }

            step *= rho;
        }

        // Line search failed
        Ok((0.0, current_energy))
    }

    /// The objective being minimized: receptor interaction plus the ligand's
    /// own internal energy.
    ///
    /// The internal term is what stops a torsion from being rotated into a
    /// self-overlapping conformation for free.
    fn calculate_energy(
        &self,
        ligand: &Molecule,
        receptor: &Molecule,
        forcefield: &dyn ForceField,
        intra_pairs: &[(usize, usize)],
    ) -> Result<f64, OptimizationError> {
        Ok(intermolecular_energy(ligand, receptor, forcefield)?
            + intramolecular_energy(ligand, intra_pairs, forcefield)?)
    }
}
