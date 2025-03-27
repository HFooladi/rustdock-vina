//! Atom representation and related functionality

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Represents different atom types supported in the forcefield
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AtomType {
    // Non-hydrogen types
    Carbon,       // C
    Nitrogen,     // N
    NitrogenH,    // NA (hydrogen bond acceptor)
    Oxygen,       // O
    OxygenH,      // OA (hydrogen bond acceptor)
    Sulfur,       // S
    SulfurH,      // SA (hydrogen bond acceptor)
    Phosphorus,   // P
    Fluorine,     // F
    Chlorine,     // Cl
    Bromine,      // Br
    Iodine,       // I
    
    // Hydrogen types
    Hydrogen,     // H
    HydrogenD,    // HD (hydrogen bond donor)
    
    // Metal types
    Zinc,         // Zn
    Calcium,      // Ca
    Manganese,    // Mn
    Magnesium,    // Mg
    Iron,         // Fe
    
    // Zinc pseudo-atom
    ZincPseudo,   // TZ
    
    // For atoms that don't match any of the above
    Unknown,
}

impl AtomType {
    /// Returns the radius of the atom type in Angstroms
    pub fn radius(&self) -> f64 {
        match self {
            AtomType::Carbon => 2.0,
            AtomType::Nitrogen => 1.75,
            AtomType::NitrogenH => 1.75,
            AtomType::Oxygen => 1.6,
            AtomType::OxygenH => 1.6,
            AtomType::Sulfur => 2.0,
            AtomType::SulfurH => 2.0,
            AtomType::Phosphorus => 2.1,
            AtomType::Fluorine => 1.54,
            AtomType::Chlorine => 2.04,
            AtomType::Bromine => 2.165,
            AtomType::Iodine => 2.36,
            AtomType::Hydrogen => 1.0,
            AtomType::HydrogenD => 1.0,
            AtomType::Zinc => 1.48,
            AtomType::Calcium => 1.98,
            AtomType::Manganese => 1.3,
            AtomType::Magnesium => 1.3,
            AtomType::Iron => 1.3,
            AtomType::ZincPseudo => 0.25, // Small radius for pseudo atom
            AtomType::Unknown => 2.0,     // Default radius
        }
    }
    
    /// Parse atom type from string representation in PDBQT format
    pub fn from_pdbqt_string(s: &str) -> Self {
        match s.trim().to_uppercase().as_str() {
            "C" => AtomType::Carbon,
            "N" => AtomType::Nitrogen,
            "NA" => AtomType::NitrogenH,
            "O" => AtomType::Oxygen,
            "OA" => AtomType::OxygenH,
            "S" => AtomType::Sulfur,
            "SA" => AtomType::SulfurH,
            "P" => AtomType::Phosphorus,
            "F" => AtomType::Fluorine,
            "CL" => AtomType::Chlorine,
            "BR" => AtomType::Bromine,
            "I" => AtomType::Iodine,
            "H" => AtomType::Hydrogen,
            "HD" => AtomType::HydrogenD,
            "ZN" => AtomType::Zinc,
            "CA" => AtomType::Calcium,
            "MN" => AtomType::Manganese,
            "MG" => AtomType::Magnesium,
            "FE" => AtomType::Iron,
            "TZ" => AtomType::ZincPseudo,
            _ => AtomType::Unknown,
        }
    }
    
    /// Convert atom type to PDBQT string
    pub fn to_pdbqt_string(&self) -> &'static str {
        match self {
            AtomType::Carbon => "C",
            AtomType::Nitrogen => "N",
            AtomType::NitrogenH => "NA",
            AtomType::Oxygen => "O",
            AtomType::OxygenH => "OA",
            AtomType::Sulfur => "S",
            AtomType::SulfurH => "SA",
            AtomType::Phosphorus => "P",
            AtomType::Fluorine => "F",
            AtomType::Chlorine => "Cl",
            AtomType::Bromine => "Br",
            AtomType::Iodine => "I",
            AtomType::Hydrogen => "H",
            AtomType::HydrogenD => "HD",
            AtomType::Zinc => "Zn",
            AtomType::Calcium => "Ca",
            AtomType::Manganese => "Mn",
            AtomType::Magnesium => "Mg",
            AtomType::Iron => "Fe",
            AtomType::ZincPseudo => "TZ",
            AtomType::Unknown => "X",
        }
    }
}

/// Represents an atom in 3D space
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Atom {
    /// Atom type 
    pub atom_type: AtomType,
    
    /// 3D coordinates (in Angstroms)
    pub coordinates: Vector3<f64>,
    
    /// Atom name from PDB format (e.g., "CA", "N", "O")
    pub name: String,
    
    /// Atom serial number from PDB
    pub serial: u32,
    
    /// Residue name this atom belongs to
    pub residue_name: String,
    
    /// Residue number this atom belongs to
    pub residue_num: u32,
    
    /// Chain identifier
    pub chain_id: char,
    
    /// Partial charge
    pub charge: f64,
    
    /// Is this atom part of a flexible side chain?
    pub is_flexible: bool,
}

impl Atom {
    /// Create a new atom
    pub fn new(
        atom_type: AtomType,
        coordinates: Vector3<f64>,
        name: String,
        serial: u32,
        residue_name: String,
        residue_num: u32,
        chain_id: char,
        charge: f64,
    ) -> Self {
        Self {
            atom_type,
            coordinates,
            name,
            serial,
            residue_name,
            residue_num,
            chain_id,
            charge,
            is_flexible: false,
        }
    }
    
    /// Calculate distance to another atom
    pub fn distance(&self, other: &Atom) -> f64 {
        (self.coordinates - other.coordinates).norm()
    }
    
    /// Check if this atom can form hydrogen bonds
    pub fn is_h_bond_donor(&self) -> bool {
        matches!(self.atom_type, AtomType::HydrogenD)
    }
    
    /// Check if this atom can accept hydrogen bonds
    pub fn is_h_bond_acceptor(&self) -> bool {
        matches!(
            self.atom_type,
            AtomType::NitrogenH | AtomType::OxygenH | AtomType::SulfurH
        )
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}({}, {}, {}) [{}]",
            self.atom_type.to_pdbqt_string(),
            self.coordinates.x,
            self.coordinates.y,
            self.coordinates.z,
            self.charge
        )
    }
}