//! Local optimization module using gradient-based methods
//!
//! This module implements BFGS-style local optimization for refining docking poses.

use nalgebra::{DVector, Matrix3, Unit, UnitQuaternion, Vector3};

use crate::forcefield::ForceField;
use crate::molecule::Molecule;
use crate::optimization::monte_carlo::{
    intramolecular_energy, intramolecular_energy_and_gradients, InterModel,
};
use crate::optimization::OptimizationError;

/// Parameters for local optimization
#[derive(Debug, Clone)]
pub struct LocalOptimizerParams {
    /// Step size for finite-difference gradients.
    ///
    /// The optimizer itself uses analytic gradients; this is only used to
    /// cross-check them.
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
        let inter = InterModel::exact(receptor);
        self.minimize_with_model(ligand, inter, forcefield, box_bounds, intra_pairs)
    }

    /// The minimization proper, scoring ligand-receptor energy through
    /// `inter` (grids during search, exact pairwise otherwise).
    pub(crate) fn minimize_with_model(
        &self,
        ligand: &mut Molecule,
        inter: InterModel,
        forcefield: &dyn ForceField,
        box_bounds: Option<(&Vector3<f64>, &Vector3<f64>)>,
        intra_pairs: &[(usize, usize)],
    ) -> Result<f64, OptimizationError> {
        let base = ligand.clone();
        let chart = Chart::new(&base);
        let n_dof = 6 + base.torsions.len();

        let mut x = DVector::zeros(n_dof);
        let mut energy = self.calculate_energy(&base, inter, forcefield, intra_pairs)?;

        // L-BFGS history
        let m = 10; // Number of corrections to store
        let mut s_history: Vec<DVector<f64>> = Vec::with_capacity(m);
        let mut y_history: Vec<DVector<f64>> = Vec::with_capacity(m);
        let mut rho_history: Vec<f64> = Vec::with_capacity(m);

        let mut grad = self.compute_gradient(&chart, inter, forcefield, &x, intra_pairs)?;

        for _ in 0..self.params.max_iterations {
            if grad.norm() < self.params.gradient_tolerance {
                break;
            }

            let direction = self.lbfgs_direction(&grad, &s_history, &y_history, &rho_history);

            let (step_size, new_energy) = self.line_search(
                &chart,
                inter,
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
            let new_grad = self.compute_gradient(&chart, inter, forcefield, &new_x, intra_pairs)?;

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

        *ligand = chart.pose(&x);
        Ok(energy)
    }

    /// Energy gradient with respect to the state, computed analytically.
    ///
    /// One pass gives every atom's `dE/dr_i`; the chain rule through
    /// [`Chart::pose`] then turns those into the translation, rotation and
    /// torsion components. This replaces `6 + n_torsions` extra energy
    /// evaluations per gradient with none.
    fn compute_gradient(
        &self,
        chart: &Chart,
        inter: InterModel,
        forcefield: &dyn ForceField,
        state: &DVector<f64>,
        intra_pairs: &[(usize, usize)],
    ) -> Result<DVector<f64>, OptimizationError> {
        let pose = chart.pose(state);
        let mut atom_grads = vec![Vector3::zeros(); pose.atoms.len()];
        inter.energy_and_gradients(&pose, forcefield, &mut atom_grads)?;
        intramolecular_energy_and_gradients(&pose, intra_pairs, forcefield, &mut atom_grads)?;
        Ok(chart.state_gradient(&pose, state, &atom_grads))
    }

    /// Forward finite-difference gradient in state space, kept to validate
    /// [`compute_gradient`](Self::compute_gradient).
    #[cfg(test)]
    fn finite_difference_gradient(
        &self,
        chart: &Chart,
        inter: InterModel,
        forcefield: &dyn ForceField,
        state: &DVector<f64>,
        intra_pairs: &[(usize, usize)],
    ) -> Result<DVector<f64>, OptimizationError> {
        let h = self.params.gradient_step;
        let mut grad = DVector::zeros(state.len());
        for i in 0..state.len() {
            let mut plus = state.clone();
            let mut minus = state.clone();
            plus[i] += h;
            minus[i] -= h;
            let e_plus =
                self.calculate_energy(&chart.pose(&plus), inter, forcefield, intra_pairs)?;
            let e_minus =
                self.calculate_energy(&chart.pose(&minus), inter, forcefield, intra_pairs)?;
            grad[i] = (e_plus - e_minus) / (2.0 * h);
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
        chart: &Chart,
        inter: InterModel,
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
            let candidate = chart.pose(&new_x);

            // Shrink the step rather than accepting a pose that leaves the box.
            let in_box = box_bounds
                .map(|(center, size)| candidate.center_within_box(center, size))
                .unwrap_or(true);

            if in_box {
                let new_energy =
                    self.calculate_energy(&candidate, inter, forcefield, intra_pairs)?;

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
        inter: InterModel,
        forcefield: &dyn ForceField,
        intra_pairs: &[(usize, usize)],
    ) -> Result<f64, OptimizationError> {
        Ok(inter.energy(ligand, forcefield)?
            + intramolecular_energy(ligand, intra_pairs, forcefield)?)
    }
}

/// The coordinate chart the local optimizer works in.
///
/// State `x = [t, phi, dtheta]` maps to positions
/// `r_i = pivot + R(phi) (p_i(dtheta) - pivot) + t`, where `p_i` applies the
/// torsion changes to the entry pose and `pivot` is the entry pose's
/// centroid. The pivot is fixed for the whole minimization: rotating about
/// the *current* centroid instead would make it depend on the torsions and
/// couple the rotation and torsion derivatives.
struct Chart<'m> {
    base: &'m Molecule,
    pivot: Vector3<f64>,
}

impl<'m> Chart<'m> {
    fn new(base: &'m Molecule) -> Self {
        let pivot = base.center().unwrap_or_else(|_| Vector3::zeros());
        Self { base, pivot }
    }

    /// Build the conformation reached by applying state `x` to the entry pose.
    ///
    /// Torsions are applied first so that the rigid-body rotation and
    /// translation act on the final internal geometry.
    fn pose(&self, x: &DVector<f64>) -> Molecule {
        let mut mol = self.base.clone();

        for i in 0..self.base.torsions.len() {
            let delta = x[6 + i];
            if delta != 0.0 {
                // A torsion that cannot be rotated (malformed tree) is skipped
                // rather than aborting the whole minimization.
                let _ = mol.rotate_torsion(i, delta);
            }
        }

        let rotation = rotation_of(x);
        let shift = self.pivot + Vector3::new(x[0], x[1], x[2]);
        for atom in &mut mol.atoms {
            atom.coordinates = rotation.transform_vector(&(atom.coordinates - self.pivot)) + shift;
        }
        mol
    }

    /// Chain rule from per-atom gradients `atom_grads` at `pose` (which must
    /// be `self.pose(x)`) to the state gradient.
    fn state_gradient(
        &self,
        pose: &Molecule,
        x: &DVector<f64>,
        atom_grads: &[Vector3<f64>],
    ) -> DVector<f64> {
        let mut grad = DVector::zeros(x.len());
        let shift = self.pivot + Vector3::new(x[0], x[1], x[2]);

        // Translation moves every atom equally; rotation contributes the
        // torque about the (shifted) pivot, mapped back through the SO(3)
        // left Jacobian because phi is a rotation vector, not an increment.
        let mut force = Vector3::zeros();
        let mut torque = Vector3::zeros();
        for (atom, g) in pose.atoms.iter().zip(atom_grads) {
            force += g;
            torque += (atom.coordinates - shift).cross(g);
        }
        let phi = Vector3::new(x[3], x[4], x[5]);
        let d_phi = left_jacobian(&phi).transpose() * torque;
        for k in 0..3 {
            grad[k] = force[k];
            grad[3 + k] = d_phi[k];
        }

        // A torsion swings its moving atoms about its bond. With a nested
        // torsion tree that motion, read in the final frame, is
        // axis x (r_i - origin) for each moving atom.
        for (k, torsion) in self.base.torsions.iter().enumerate() {
            let Some(bond) = self.base.bonds.get(torsion.bond_idx) else {
                continue;
            };
            let (Some(a1), Some(a2)) = (
                pose.atoms.get(bond.atom1_idx),
                pose.atoms.get(bond.atom2_idx),
            ) else {
                continue;
            };
            let axis = a2.coordinates - a1.coordinates;
            if axis.norm() < 1e-10 || torsion.moving_atoms.iter().any(|&i| i >= pose.atoms.len()) {
                continue;
            }
            let axis = axis.normalize();
            grad[6 + k] = torsion
                .moving_atoms
                .iter()
                .map(|&i| {
                    atom_grads[i].dot(&axis.cross(&(pose.atoms[i].coordinates - a1.coordinates)))
                })
                .sum();
        }

        grad
    }
}

/// The rotation encoded by the rotation-vector part of the state.
fn rotation_of(x: &DVector<f64>) -> UnitQuaternion<f64> {
    let phi = Vector3::new(x[3], x[4], x[5]);
    let angle = phi.norm();
    if angle > 1e-10 {
        UnitQuaternion::from_axis_angle(&Unit::new_normalize(phi), angle)
    } else {
        UnitQuaternion::identity()
    }
}

/// Left Jacobian of SO(3) at rotation vector `phi`: how a change in `phi`
/// shows up as an infinitesimal rotation applied on the left.
fn left_jacobian(phi: &Vector3<f64>) -> Matrix3<f64> {
    let theta = phi.norm();
    let k = phi.cross_matrix();
    let (a, b) = if theta < 1e-4 {
        // Series expansions of the coefficients below.
        (
            0.5 - theta * theta / 24.0,
            1.0 / 6.0 - theta * theta / 120.0,
        )
    } else {
        let t2 = theta * theta;
        (
            (1.0 - theta.cos()) / t2,
            (theta - theta.sin()) / (t2 * theta),
        )
    };
    Matrix3::identity() + k * a + k * k * b
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::vina::VinaForceField;
    use crate::grid::GridSet;
    use crate::io::parse_pdbqt;
    use rand::prelude::*;
    use rand::rngs::StdRng;
    use std::path::PathBuf;

    fn benchmark(name: &str) -> Molecule {
        parse_pdbqt(
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("benchmark_data")
                .join(name),
        )
        .unwrap()
    }

    /// The analytic state gradient must match finite differences on a real
    /// flexible ligand, for both grid and exact scoring, at states with
    /// non-trivial rotation and torsion offsets. This is also the check that
    /// the parsed torsion tree is nested, which the torsion formula assumes.
    #[test]
    fn test_analytic_gradient_matches_finite_difference() {
        let receptor = benchmark("receptor.pdbqt");
        let ligand = benchmark("ligand.pdbqt");
        let ff = VinaForceField::default();
        let pairs = ligand.intramolecular_pairs();
        let center = ligand.center().unwrap();
        let grids = GridSet::build(
            &ligand,
            &receptor,
            &ff,
            &center,
            &Vector3::repeat(16.0),
            4.0,
            0.375,
        )
        .unwrap()
        .unwrap();

        let optimizer = LocalOptimizer::with_params(LocalOptimizerParams {
            gradient_step: 1e-5,
            ..LocalOptimizerParams::default()
        });
        let chart = Chart::new(&ligand);
        let n_dof = 6 + ligand.torsions.len();
        let mut rng = StdRng::seed_from_u64(3);

        for model in [
            InterModel::exact(&receptor),
            InterModel::with_grids(&receptor, &grids),
        ] {
            for _ in 0..4 {
                let x = DVector::from_fn(n_dof, |i, _| {
                    let scale = if i < 3 { 0.5 } else { 0.6 };
                    (rng.gen::<f64>() - 0.5) * 2.0 * scale
                });
                let analytic = optimizer
                    .compute_gradient(&chart, model, &ff, &x, &pairs)
                    .unwrap();
                let numeric = optimizer
                    .finite_difference_gradient(&chart, model, &ff, &x, &pairs)
                    .unwrap();
                let err = (&analytic - &numeric).norm() / numeric.norm().max(1.0);
                assert!(
                    err < 1e-3,
                    "relative error {err}\nanalytic {analytic}\nnumeric {numeric}"
                );
            }
        }
    }

    #[test]
    fn test_left_jacobian_matches_finite_difference() {
        let phi = Vector3::new(0.7, -0.4, 1.1);
        let v = Vector3::new(0.3, 2.0, -1.0);
        let r = |p: &Vector3<f64>| {
            let x = DVector::from_vec(vec![0.0, 0.0, 0.0, p.x, p.y, p.z]);
            rotation_of(&x).transform_vector(&v)
        };
        let jac = left_jacobian(&phi);
        let rv = r(&phi);
        let h = 1e-6;
        for k in 0..3 {
            let mut dp = Vector3::zeros();
            dp[k] = h;
            let numeric = (r(&(phi + dp)) - r(&(phi - dp))) / (2.0 * h);
            let analytic = (jac * dp / h).cross(&rv);
            assert!((numeric - analytic).norm() < 1e-6, "column {k}");
        }
    }
}
