//! Implementation of the AutoDock4 scoring function

use crate::atom::{Atom, AtomType};
use crate::forcefield::{ForceField, ForceFieldError};
use nalgebra::Vector3;
use std::collections::HashMap;

/// Parameters for AutoDock4 forcefield
#[derive(Debug, Clone)]
pub struct AD4Params {
    // Weights for each component of the scoring function
    pub weight_vdw: f64,
    pub weight_hbond: f64,
    pub weight_elec: f64,
    pub weight_desolv: f64,
    pub weight_tors: f64,

    // vdW parameters
    pub vdw_coeffs: HashMap<(AtomType, AtomType), (f64, f64)>, // (Rii, epsii)

    // Hydrogen bond parameters
    pub hbond_coeffs: HashMap<(AtomType, AtomType), (f64, f64)>, // (Rij_hb, epsij_hb)

    // Desolvation parameters
    pub atom_volumes: HashMap<AtomType, f64>,
    pub atom_solvpars: HashMap<AtomType, f64>,
}

impl Default for AD4Params {
    fn default() -> Self {
        let mut vdw_coeffs = HashMap::new();
        let mut hbond_coeffs = HashMap::new();
        let mut atom_volumes = HashMap::new();
        let mut atom_solvpars = HashMap::new();

        // Initialize with some default values
        // In a real implementation, these would be loaded from a parameter file

        // Carbon
        atom_volumes.insert(AtomType::Carbon, 33.5103);
        atom_solvpars.insert(AtomType::Carbon, -0.00143);

        // Nitrogen
        atom_volumes.insert(AtomType::Nitrogen, 22.4493);
        atom_volumes.insert(AtomType::NitrogenH, 22.4493);
        atom_solvpars.insert(AtomType::Nitrogen, -0.00162);
        atom_solvpars.insert(AtomType::NitrogenH, -0.00162);

        // Oxygen
        atom_volumes.insert(AtomType::Oxygen, 17.1573);
        atom_volumes.insert(AtomType::OxygenH, 17.1573);
        atom_solvpars.insert(AtomType::Oxygen, -0.00251);
        atom_solvpars.insert(AtomType::OxygenH, -0.00251);

        // Hydrogen
        atom_volumes.insert(AtomType::Hydrogen, 0.0);
        atom_volumes.insert(AtomType::HydrogenD, 0.0);
        atom_solvpars.insert(AtomType::Hydrogen, 0.00051);
        atom_solvpars.insert(AtomType::HydrogenD, 0.00051);

        // Example vdW parameters (carbon-carbon)
        vdw_coeffs.insert(
            (AtomType::Carbon, AtomType::Carbon),
            (4.0, 0.15), // (Rii, epsii)
        );

        // Example H-bond parameters (hydrogen donor - oxygen acceptor)
        hbond_coeffs.insert(
            (AtomType::HydrogenD, AtomType::OxygenH),
            (1.9, 5.0), // (Rij_hb, epsij_hb)
        );

        Self {
            weight_vdw: 0.1662,
            weight_hbond: 0.1209,
            weight_elec: 0.1406,
            weight_desolv: 0.1322,
            weight_tors: 0.2983,

            vdw_coeffs,
            hbond_coeffs,
            atom_volumes,
            atom_solvpars,
        }
    }
}

/// Implementation of the AutoDock4 scoring function
#[derive(Debug, Clone)]
pub struct AD4ForceField {
    pub params: AD4Params,
}

impl Default for AD4ForceField {
    fn default() -> Self {
        Self {
            params: AD4Params::default(),
        }
    }
}

impl ForceField for AD4ForceField {
    fn name(&self) -> &'static str {
        "AutoDock4"
    }

    fn atom_pair_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError> {
        // Calculate all components
        let vdw = self.vdw_energy(atom1, atom2, distance)?;
        let hbond = self.hbond_energy(atom1, atom2, distance)?;
        let elec = self.electrostatic_energy(atom1, atom2, distance)?;
        let desolv = self.desolvation_energy(atom1, atom2, distance)?;

        // Apply weights and sum
        let energy = self.params.weight_vdw * vdw
            + self.params.weight_hbond * hbond
            + self.params.weight_elec * elec
            + self.params.weight_desolv * desolv;

        Ok(energy)
    }

    fn vdw_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError> {
        // Get vdW parameters
        let type1 = atom1.atom_type;
        let type2 = atom2.atom_type;

        // Try both orders of atom types
        let params = self
            .params
            .vdw_coeffs
            .get(&(type1, type2))
            .or_else(|| self.params.vdw_coeffs.get(&(type2, type1)));

        let (r_sum, eps) = match params {
            Some(p) => *p,
            None => {
                // Calculate arithmetic mean of radii and geometric mean of well depths
                let r1 = type1.radius();
                let r2 = type2.radius();
                let eps1 = self.get_epsilon(atom1)?;
                let eps2 = self.get_epsilon(atom2)?;

                ((r1 + r2) / 2.0, (eps1 * eps2).sqrt())
            }
        };

        // Calculate vdW energy using 12-6 Lennard-Jones potential
        if distance < 0.01 {
            return Ok(1000.0); // Arbitrary high repulsion for very close atoms
        }

        let r_ratio = r_sum / distance;
        let r6 = r_ratio.powi(6);
        let r12 = r6 * r6;

        let energy = eps * (r12 - 2.0 * r6);

        Ok(energy)
    }

    fn electrostatic_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError> {
        // AutoDock4 uses a distance-dependent dielectric model
        if distance < 0.01 {
            return Ok(0.0); // Avoid division by zero
        }

        // Distance-dependent dielectric function from AutoDock4
        let epsilon = 4.0 * distance; // Simplified version

        // Coulomb's law
        let energy = 332.0 * atom1.charge * atom2.charge / (distance * epsilon);

        Ok(energy)
    }

    fn hbond_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError> {
        // Check if we have a hydrogen bond donor and acceptor
        let (donor_atom, acceptor_atom) = if atom1.is_h_bond_donor() && atom2.is_h_bond_acceptor() {
            (atom1, atom2)
        } else if atom2.is_h_bond_donor() && atom1.is_h_bond_acceptor() {
            (atom2, atom1)
        } else {
            // No hydrogen bond possible
            return Ok(0.0);
        };

        // Get parameters
        let donor_type = donor_atom.atom_type;
        let acceptor_type = acceptor_atom.atom_type;

        let params = self
            .params
            .hbond_coeffs
            .get(&(donor_type, acceptor_type))
            .or_else(|| self.params.hbond_coeffs.get(&(acceptor_type, donor_type)));

        let (optimal_dist, well_depth) = match params {
            Some(p) => *p,
            None => (1.9, 1.0), // Default values
        };

        // Calculate H-bond energy using 12-10 potential
        if distance < 0.01 {
            return Ok(0.0); // No H-bond at very short range
        }

        let r_ratio = optimal_dist / distance;
        let r10 = r_ratio.powi(10);
        let r12 = r_ratio.powi(12);

        let energy = well_depth * (5.0 * r12 - 6.0 * r10);

        Ok(energy)
    }

    fn desolvation_energy(
        &self,
        atom1: &Atom,
        atom2: &Atom,
        distance: f64,
    ) -> Result<f64, ForceFieldError> {
        // Get desolvation parameters
        let vol1 = self
            .params
            .atom_volumes
            .get(&atom1.atom_type)
            .copied()
            .unwrap_or(0.0);
        let vol2 = self
            .params
            .atom_volumes
            .get(&atom2.atom_type)
            .copied()
            .unwrap_or(0.0);
        let solv1 = self
            .params
            .atom_solvpars
            .get(&atom1.atom_type)
            .copied()
            .unwrap_or(0.0);
        let solv2 = self
            .params
            .atom_solvpars
            .get(&atom2.atom_type)
            .copied()
            .unwrap_or(0.0);

        // Calculate desolvation using AutoDock4 formula
        let sigma = 3.6; // Default value from AutoDock4
        let r_minus_sigma = distance - sigma;

        let exponent = -0.5 * r_minus_sigma * r_minus_sigma / (sigma * sigma);
        let desolv_energy = (solv1 * vol2 + solv2 * vol1) * exponent.exp();

        Ok(desolv_energy)
    }

    fn atom_charge_energy(
        &self,
        atom: &Atom,
        charge: f64,
        position: &Vector3<f64>,
    ) -> Result<f64, ForceFieldError> {
        let distance = (atom.coordinates - position).norm();

        // Early return if distance is too small
        if distance < 0.01 {
            return Ok(0.0);
        }

        // Distance-dependent dielectric
        let epsilon = 4.0 * distance; // Simplified

        // Coulomb's law
        let energy = 332.0 * atom.charge * charge / (distance * epsilon);

        Ok(energy)
    }
}

// Additional methods specific to AutoDock4
impl AD4ForceField {
    /// Calculate the torsional entropy penalty
    pub fn torsional_entropy(&self, num_rotatable_bonds: usize) -> f64 {
        self.params.weight_tors * num_rotatable_bonds as f64
    }
}
