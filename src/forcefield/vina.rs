//! Implementation of the Vina scoring function
//!
//! This implementation follows the original AutoDock Vina scoring function
//! as described in Trott & Olson, J. Comput. Chem. 2010.

use crate::atom::{Atom, AtomType};
use crate::forcefield::{ForceField, ForceFieldError};
use nalgebra::Vector3;

/// Parameters for Vina forcefield
/// Weights are from the original Vina paper (Trott & Olson, 2010)
#[derive(Debug, Clone)]
pub struct VinaParams {
    /// Weight for first Gaussian term (steric, short-range attractive)
    pub weight_gauss1: f64,
    /// Weight for second Gaussian term (steric, longer-range attractive)
    pub weight_gauss2: f64,
    /// Weight for repulsion term (steric clash penalty)
    pub weight_repulsion: f64,
    /// Weight for hydrophobic term
    pub weight_hydrophobic: f64,
    /// Weight for hydrogen bonding term
    pub weight_hydrogen: f64,
    /// Weight for rotatable bond penalty (entropy)
    pub weight_rot: f64,
}

impl Default for VinaParams {
    fn default() -> Self {
        // Exact weights from the original AutoDock Vina
        // Source: Trott & Olson, J. Comput. Chem. 2010
        Self {
            weight_gauss1: -0.035579,
            weight_gauss2: -0.005156,
            weight_repulsion: 0.840245,
            weight_hydrophobic: -0.035069,
            weight_hydrogen: -0.587439,
            weight_rot: 0.05846,
        }
    }
}

/// Implementation of the Vina scoring function
#[derive(Debug, Clone, Default)]
pub struct VinaForceField {
    pub params: VinaParams,
}

impl VinaForceField {
    /// Create a new VinaForceField with default parameters
    pub fn new() -> Self {
        Self::default()
    }

    /// Calculate the surface distance between two atoms
    /// This is distance minus the sum of van der Waals radii
    #[inline]
    fn surface_distance(atom1: &Atom, atom2: &Atom, distance: f64) -> f64 {
        let r1 = atom1.atom_type.radius();
        let r2 = atom2.atom_type.radius();
        distance - r1 - r2
    }

    /// Check if an atom type is hydrophobic (carbon-like)
    #[inline]
    fn is_hydrophobic(atom_type: AtomType) -> bool {
        matches!(
            atom_type,
            AtomType::Carbon
                | AtomType::Fluorine
                | AtomType::Chlorine
                | AtomType::Bromine
                | AtomType::Iodine
        )
    }

    /// Check if an atom type can be a hydrogen bond donor
    /// In PDBQT, NA/OA/SA are marked as acceptors but in practice they often
    /// have attached hydrogens and can also act as donors
    #[inline]
    fn is_hbond_donor(atom_type: AtomType) -> bool {
        // NitrogenH, OxygenH, SulfurH correspond to PDBQT NA, OA, SA
        // which can be both donors and acceptors (like Vina's N_DA, O_DA types)
        matches!(
            atom_type,
            AtomType::NitrogenH | AtomType::OxygenH | AtomType::SulfurH
        )
    }

    /// Check if an atom type can be a hydrogen bond acceptor
    #[inline]
    fn is_hbond_acceptor(atom_type: AtomType) -> bool {
        // All nitrogen/oxygen/sulfur types can accept H-bonds
        matches!(
            atom_type,
            AtomType::Nitrogen
                | AtomType::NitrogenH
                | AtomType::Oxygen
                | AtomType::OxygenH
                | AtomType::Sulfur
                | AtomType::SulfurH
        )
    }

    /// Calculate the conformational entropy penalty for rotatable bonds
    /// This is added to the final score: penalty = weight_rot * N_rot
    pub fn rotatable_bond_penalty(&self, num_rotatable_bonds: usize) -> f64 {
        self.params.weight_rot * num_rotatable_bonds as f64
    }

    /// Calculate the full interaction energy and apply torsional normalization
    /// affinity = (intermolecular_energy) / (1 + weight_rot * N_rot)
    pub fn calculate_affinity(
        &self,
        intermolecular_energy: f64,
        num_rotatable_bonds: usize,
    ) -> f64 {
        let divisor = 1.0 + self.params.weight_rot * num_rotatable_bonds as f64;
        intermolecular_energy / divisor
    }
}

impl ForceField for VinaForceField {
    fn name(&self) -> &'static str {
        "Vina"
    }

    fn atom_pair_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError> {
        // Cutoff distance (8 Angstroms)
        if !(0.01..=8.0).contains(&distance) {
            return Ok(0.0);
        }

        // Calculate surface distance (d = r - r1 - r2)
        let d = Self::surface_distance(atom1, atom2, distance);

        // Calculate all interaction components
        let gauss1 = self.gauss1_term(d);
        let gauss2 = self.gauss2_term(d);
        let repulsion = self.repulsion_term(d);
        let hydrophobic = self.hydrophobic_term(atom1, atom2, d);
        let hbond = self.hbond_term(atom1, atom2, d);

        Ok(gauss1 + gauss2 + repulsion + hydrophobic + hbond)
    }

    fn vdw_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError> {
        if !(0.01..=8.0).contains(&distance) {
            return Ok(0.0);
        }

        let d = Self::surface_distance(atom1, atom2, distance);

        // van der Waals in Vina is gauss1 + gauss2 + repulsion
        let gauss1 = self.gauss1_term(d);
        let gauss2 = self.gauss2_term(d);
        let repulsion = self.repulsion_term(d);

        Ok(gauss1 + gauss2 + repulsion)
    }

    fn electrostatic_energy(
        &self,
        _atom1: &Atom,
        _atom2: &Atom,
        _distance: f64,
    ) -> Result<f64, ForceFieldError> {
        // Original Vina does not have an explicit electrostatic term
        // Electrostatics are implicitly captured in the empirical weights
        Ok(0.0)
    }

    fn hbond_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError> {
        if !(0.01..=8.0).contains(&distance) {
            return Ok(0.0);
        }

        let d = Self::surface_distance(atom1, atom2, distance);
        Ok(self.hbond_term(atom1, atom2, d))
    }

    fn hydrophobic_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError> {
        if !(0.01..=8.0).contains(&distance) {
            return Ok(0.0);
        }

        let d = Self::surface_distance(atom1, atom2, distance);
        Ok(self.hydrophobic_term(atom1, atom2, d))
    }

    fn desolvation_energy(
        &self,
        _atom1: &Atom,
        _atom2: &Atom,
        _distance: f64,
    ) -> Result<f64, ForceFieldError> {
        // Original Vina does not have an explicit desolvation term
        // Desolvation is implicitly captured in the empirical weights
        Ok(0.0)
    }

    fn atom_charge_energy(
        &self,
        _atom: &Atom,
        _charge: f64,
        _position: &Vector3<f64>,
    ) -> Result<f64, ForceFieldError> {
        // Not used in Vina scoring
        Ok(0.0)
    }

    /// Apply Vina's torsional normalization: affinity = energy / (1 + weight_rot * N_rot)
    fn apply_torsional_normalization(
        &self,
        intermolecular_energy: f64,
        num_rotatable_bonds: usize,
    ) -> f64 {
        self.calculate_affinity(intermolecular_energy, num_rotatable_bonds)
    }
}

// Private implementation of individual scoring terms
impl VinaForceField {
    /// First Gaussian term: attractive, short-range
    /// gauss1 = w1 * exp(-(d/0.5)^2) where d is surface distance
    #[inline]
    fn gauss1_term(&self, d: f64) -> f64 {
        // Vina formula: exp(-(d/0.5)^2) = exp(-4*d^2)
        self.params.weight_gauss1 * (-(d / 0.5).powi(2)).exp()
    }

    /// Second Gaussian term: attractive, longer-range
    /// gauss2 = w2 * exp(-((d-3)/2)^2)
    #[inline]
    fn gauss2_term(&self, d: f64) -> f64 {
        // Vina formula: exp(-((d-3)/2)^2)
        self.params.weight_gauss2 * (-((d - 3.0) / 2.0).powi(2)).exp()
    }

    /// Repulsion term: penalty for steric clash
    /// repulsion = w3 * d^2 if d < 0, else 0
    #[inline]
    fn repulsion_term(&self, d: f64) -> f64 {
        if d < 0.0 {
            self.params.weight_repulsion * d.powi(2)
        } else {
            0.0
        }
    }

    /// Hydrophobic term: attractive for hydrophobic-hydrophobic contacts
    /// Uses piecewise linear function:
    ///   if d < 0.5: 1
    ///   if 0.5 <= d < 1.5: linear interpolation to 0
    ///   if d >= 1.5: 0
    #[inline]
    fn hydrophobic_term(&self, atom1: &Atom, atom2: &Atom, d: f64) -> f64 {
        // Both atoms must be hydrophobic
        if !Self::is_hydrophobic(atom1.atom_type) || !Self::is_hydrophobic(atom2.atom_type) {
            return 0.0;
        }

        let h = if d < 0.5 {
            1.0
        } else if d < 1.5 {
            1.5 - d
        } else {
            0.0
        };

        self.params.weight_hydrophobic * h
    }

    /// Hydrogen bond term: attractive for donor-acceptor pairs
    /// Uses piecewise linear function:
    ///   if d < -0.7: 1
    ///   if -0.7 <= d < 0: linear interpolation
    ///   if d >= 0: 0
    #[inline]
    fn hbond_term(&self, atom1: &Atom, atom2: &Atom, d: f64) -> f64 {
        // Check for donor-acceptor pair
        let is_hbond_pair = (Self::is_hbond_donor(atom1.atom_type)
            && Self::is_hbond_acceptor(atom2.atom_type))
            || (Self::is_hbond_donor(atom2.atom_type) && Self::is_hbond_acceptor(atom1.atom_type));

        if !is_hbond_pair {
            return 0.0;
        }

        let h = if d < -0.7 {
            1.0
        } else if d < 0.0 {
            -d / 0.7
        } else {
            0.0
        };

        self.params.weight_hydrogen * h
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    fn create_test_atom(atom_type: AtomType, x: f64, y: f64, z: f64) -> Atom {
        Atom::new(
            atom_type,
            Vector3::new(x, y, z),
            "TEST".to_string(),
            1,                 // serial
            "TST".to_string(), // residue_name
            1,                 // residue_num
            'A',               // chain_id
            0.0,               // charge
        )
    }

    #[test]
    fn test_gauss1_at_surface() {
        let ff = VinaForceField::new();
        // At surface distance d=0, gauss1 should be at maximum
        let energy = ff.gauss1_term(0.0);
        assert!(energy < 0.0, "Gauss1 should be negative (attractive)");
        assert!((energy - ff.params.weight_gauss1).abs() < 1e-6);
    }

    #[test]
    fn test_repulsion_only_when_overlapping() {
        let ff = VinaForceField::new();

        // No repulsion when d >= 0
        assert_eq!(ff.repulsion_term(0.0), 0.0);
        assert_eq!(ff.repulsion_term(1.0), 0.0);

        // Repulsion when d < 0 (atoms overlapping)
        let rep = ff.repulsion_term(-0.5);
        assert!(rep > 0.0, "Repulsion should be positive");
    }

    #[test]
    fn test_hydrophobic_term() {
        let ff = VinaForceField::new();
        let c1 = create_test_atom(AtomType::Carbon, 0.0, 0.0, 0.0);
        let c2 = create_test_atom(AtomType::Carbon, 4.0, 0.0, 0.0);

        // Carbon-carbon should have hydrophobic interaction
        let energy = ff.hydrophobic_term(&c1, &c2, 0.3);
        assert!(energy < 0.0, "Hydrophobic should be attractive (negative)");
    }

    #[test]
    fn test_hbond_term() {
        let ff = VinaForceField::new();
        let donor = create_test_atom(AtomType::NitrogenH, 0.0, 0.0, 0.0);
        let acceptor = create_test_atom(AtomType::Oxygen, 2.8, 0.0, 0.0);

        // Should have H-bond interaction at appropriate distance
        let energy = ff.hbond_term(&donor, &acceptor, -0.5);
        assert!(energy < 0.0, "H-bond should be attractive (negative)");
    }

    #[test]
    fn test_torsional_penalty() {
        let ff = VinaForceField::new();

        // 6 rotatable bonds (like imatinib)
        let penalty = ff.rotatable_bond_penalty(6);
        assert!((penalty - 0.35076).abs() < 0.001);
    }
}
