//! The Vinardo scoring function.
//!
//! Vinardo (Quiroga & Villarreal, PLoS ONE 2016) keeps Vina's form and atom
//! typing but drops the second Gaussian, changes the radii and reweights the
//! rest. This follows AutoDock Vina 1.2.7's implementation (`potentials.h`,
//! `scoring_function.h`, `xs_vinardo_vdw_radii`).

use crate::atom::{Atom, AtomType};
use crate::forcefield::{ForceField, ForceFieldError};
use nalgebra::Vector3;

/// Interaction cutoff in Angstroms, on the centre-to-centre distance.
pub const CUTOFF: f64 = 8.0;

/// Vinardo's weights, Vina 1.2.7's defaults.
#[derive(Debug, Clone)]
pub struct VinardoParams {
    /// Gaussian attraction (width 0.8)
    pub weight_gauss: f64,
    /// Steric repulsion
    pub weight_repulsion: f64,
    /// Hydrophobic contact
    pub weight_hydrophobic: f64,
    /// Hydrogen bond
    pub weight_hydrogen: f64,
    /// Rotatable-bond weight of the torsional normalization
    pub weight_rot: f64,
}

impl Default for VinardoParams {
    fn default() -> Self {
        Self {
            weight_gauss: -0.045,
            weight_repulsion: 0.8,
            weight_hydrophobic: -0.035,
            weight_hydrogen: -0.6,
            weight_rot: 0.05846,
        }
    }
}

/// The Vinardo scoring function.
#[derive(Debug, Clone, Default)]
pub struct VinardoForceField {
    pub params: VinardoParams,
}

impl VinardoForceField {
    pub fn new() -> Self {
        Self::default()
    }

    /// Vinardo's van der Waals radius for an atom type (Vina's
    /// `xs_vinardo_vdw_radii`): larger carbon, smaller nitrogen and oxygen
    /// than Vina's.
    pub fn radius(atom_type: AtomType) -> f64 {
        match atom_type {
            AtomType::Carbon => 2.0,
            AtomType::Nitrogen | AtomType::NitrogenH => 1.7,
            AtomType::Oxygen | AtomType::OxygenH => 1.6,
            AtomType::Sulfur | AtomType::SulfurH => 2.0,
            AtomType::Phosphorus => 2.1,
            AtomType::Fluorine => 1.5,
            AtomType::Chlorine => 1.8,
            AtomType::Bromine => 2.0,
            AtomType::Iodine => 2.2,
            t if t.is_metal() => 1.2,
            t => t.radius(),
        }
    }

    /// Surface distance `r - R1 - R2` with Vinardo's radii.
    #[inline]
    fn surface_distance(atom1: &Atom, atom2: &Atom, distance: f64) -> f64 {
        distance - Self::radius(atom1.atom_type) - Self::radius(atom2.atom_type)
    }

    #[inline]
    fn scored(atom1: &Atom, atom2: &Atom, distance: f64) -> bool {
        !atom1.atom_type.is_hydrogen()
            && !atom2.atom_type.is_hydrogen()
            && (0.01..CUTOFF).contains(&distance)
    }

    #[inline]
    fn hydrophobic_pair(atom1: &Atom, atom2: &Atom) -> bool {
        atom1.is_hydrophobic && atom2.is_hydrophobic
    }

    #[inline]
    fn hbond_pair(atom1: &Atom, atom2: &Atom) -> bool {
        (atom1.is_hbond_donor && atom2.atom_type.is_hbond_acceptor())
            || (atom2.is_hbond_donor && atom1.atom_type.is_hbond_acceptor())
    }

    #[inline]
    fn gauss(&self, d: f64) -> f64 {
        self.params.weight_gauss * (-(d / 0.8).powi(2)).exp()
    }

    #[inline]
    fn repulsion(&self, d: f64) -> f64 {
        if d < 0.0 {
            self.params.weight_repulsion * d * d
        } else {
            0.0
        }
    }

    /// 1 for `d <= 0`, falling linearly to 0 at `d = 2.5`.
    #[inline]
    fn hydrophobic(&self, d: f64) -> f64 {
        self.params.weight_hydrophobic * ((2.5 - d) / 2.5).clamp(0.0, 1.0)
    }

    /// 1 for `d <= -0.6`, falling linearly to 0 at `d = 0`.
    #[inline]
    fn hbond(&self, d: f64) -> f64 {
        self.params.weight_hydrogen * (-d / 0.6).clamp(0.0, 1.0)
    }
}

impl ForceField for VinardoForceField {
    fn name(&self) -> &'static str {
        "Vinardo"
    }

    fn atom_pair_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError> {
        if !Self::scored(atom1, atom2, distance) {
            return Ok(0.0);
        }
        let d = Self::surface_distance(atom1, atom2, distance);
        let mut e = self.gauss(d) + self.repulsion(d);
        if Self::hydrophobic_pair(atom1, atom2) {
            e += self.hydrophobic(d);
        }
        if Self::hbond_pair(atom1, atom2) {
            e += self.hbond(d);
        }
        Ok(e)
    }

    fn atom_pair_derivative(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError> {
        if !Self::scored(atom1, atom2, distance) {
            return Ok(0.0);
        }
        // Every term is a function of d = r - R1 - R2, so dE/dr = dE/dd.
        let d = Self::surface_distance(atom1, atom2, distance);
        let mut de = self.gauss(d) * (-2.0 * d / (0.8 * 0.8));
        if d < 0.0 {
            de += 2.0 * self.params.weight_repulsion * d;
        }
        if Self::hydrophobic_pair(atom1, atom2) && (0.0..2.5).contains(&d) {
            de -= self.params.weight_hydrophobic / 2.5;
        }
        if Self::hbond_pair(atom1, atom2) && (-0.6..0.0).contains(&d) {
            de -= self.params.weight_hydrogen / 0.6;
        }
        Ok(de)
    }

    fn vdw_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError> {
        if !Self::scored(atom1, atom2, distance) {
            return Ok(0.0);
        }
        let d = Self::surface_distance(atom1, atom2, distance);
        Ok(self.gauss(d) + self.repulsion(d))
    }

    fn electrostatic_energy(
        &self,
        _atom1: &Atom,
        _atom2: &Atom,
        _distance: f64,
    ) -> Result<f64, ForceFieldError> {
        // Like Vina, Vinardo has no explicit electrostatics.
        Ok(0.0)
    }

    fn hbond_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError> {
        if !Self::scored(atom1, atom2, distance) || !Self::hbond_pair(atom1, atom2) {
            return Ok(0.0);
        }
        Ok(self.hbond(Self::surface_distance(atom1, atom2, distance)))
    }

    fn desolvation_energy(
        &self,
        _atom1: &Atom,
        _atom2: &Atom,
        _distance: f64,
    ) -> Result<f64, ForceFieldError> {
        Ok(0.0)
    }

    fn hydrophobic_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError> {
        if !Self::scored(atom1, atom2, distance) || !Self::hydrophobic_pair(atom1, atom2) {
            return Ok(0.0);
        }
        Ok(self.hydrophobic(Self::surface_distance(atom1, atom2, distance)))
    }

    fn atom_charge_energy(
        &self,
        _atom: &Atom,
        _charge: f64,
        _position: &Vector3<f64>,
    ) -> Result<f64, ForceFieldError> {
        Ok(0.0)
    }

    /// Vina's torsional normalization, which Vinardo shares:
    /// affinity = energy / (1 + weight_rot * num_tors).
    fn apply_torsional_normalization(&self, intermolecular_energy: f64, num_tors: f64) -> f64 {
        intermolecular_energy / (1.0 + self.params.weight_rot * num_tors)
    }

    /// Like Vina's, Vinardo's terms depend only on `xs` type and distance.
    fn grid_cutoff(&self) -> Option<f64> {
        Some(CUTOFF)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom(atom_type: AtomType, hydrophobic: bool, donor: bool) -> Atom {
        let mut a = Atom::new(
            atom_type,
            Vector3::zeros(),
            "X".to_string(),
            1,
            "LIG".to_string(),
            1,
            'A',
            0.0,
        );
        a.is_hydrophobic = hydrophobic;
        a.is_hbond_donor = donor;
        a
    }

    #[test]
    fn test_terms_at_reference_points() {
        let ff = VinardoForceField::new();
        let c = atom(AtomType::Carbon, true, false);
        // C-C: R1 + R2 = 4.0, so r = 4.0 is touching (d = 0): full gauss
        // and full hydrophobic, no repulsion.
        let e = ff.atom_pair_energy(&c, &c, 4.0).unwrap();
        assert!((e - (-0.045 - 0.035)).abs() < 1e-12);
        // d = 2.5: hydrophobic has just reached zero.
        let e = ff.atom_pair_energy(&c, &c, 6.5).unwrap();
        assert!((e - -0.045 * (-(2.5f64 / 0.8).powi(2)).exp()).abs() < 1e-12);
        // N donor to O acceptor at d = -0.6 (r = 3.3 - 0.6): full H-bond.
        let n = atom(AtomType::Nitrogen, false, true);
        let o = atom(AtomType::OxygenH, false, false);
        let d: f64 = -0.6;
        let expected = -0.045 * (-(d / 0.8).powi(2)).exp() + 0.8 * d * d - 0.6;
        let e = ff.atom_pair_energy(&n, &o, 3.3 + d).unwrap();
        assert!((e - expected).abs() < 1e-12);
    }

    #[test]
    fn test_derivative_matches_finite_difference() {
        let ff = VinardoForceField::new();
        let pairs = [
            (
                atom(AtomType::Carbon, true, false),
                atom(AtomType::Carbon, true, false),
            ),
            (
                atom(AtomType::Nitrogen, false, true),
                atom(AtomType::OxygenH, false, false),
            ),
            (
                atom(AtomType::Chlorine, true, false),
                atom(AtomType::Carbon, false, false),
            ),
        ];
        for (a, b) in &pairs {
            let mut r = 1.5;
            while r < 7.9 {
                let h = 1e-6;
                let numeric = (ff.atom_pair_energy(a, b, r + h).unwrap()
                    - ff.atom_pair_energy(a, b, r - h).unwrap())
                    / (2.0 * h);
                let analytic = ff.atom_pair_derivative(a, b, r).unwrap();
                assert!(
                    (numeric - analytic).abs() < 1e-5,
                    "{:?}-{:?} at r={r}: {analytic} vs {numeric}",
                    a.atom_type,
                    b.atom_type
                );
                r += 0.0137;
            }
        }
    }
}
