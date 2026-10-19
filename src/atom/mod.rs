//! Atom representation and related functionality

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Represents different atom types supported in the forcefield
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AtomType {
    // Non-hydrogen types
    Carbon,     // C
    Nitrogen,   // N
    NitrogenH,  // NA (hydrogen bond acceptor)
    Oxygen,     // O
    OxygenH,    // OA (hydrogen bond acceptor)
    Sulfur,     // S
    SulfurH,    // SA (hydrogen bond acceptor)
    Phosphorus, // P
    Fluorine,   // F
    Chlorine,   // Cl
    Bromine,    // Br
    Iodine,     // I

    // Hydrogen types
    Hydrogen,  // H
    HydrogenD, // HD (hydrogen bond donor)

    // Metal types
    Zinc,      // Zn
    Calcium,   // Ca
    Manganese, // Mn
    Magnesium, // Mg
    Iron,      // Fe

    // Metals that are not AutoDock 4 types; Vina scores them all as generic
    // metals (`Met_D`), matching the names case-sensitively.
    Sodium,    // Na
    Potassium, // K
    Copper,    // Cu
    Nickel,    // Ni
    Cobalt,    // Co
    Cadmium,   // Cd
    Mercury,   // Hg
    Uranium,   // U

    // Zinc pseudo-atom
    ZincPseudo, // TZ

    // Macrocycle closure, as Meeko writes a ring it opened for flexible
    // docking: the carbon on each side of the broken bond (`CG0`-`CG3`) and
    // the dummy standing in for its partner (`G0`-`G3`). The number pairs
    // each dummy with the carbon it is glued back onto.
    CarbonClosure(u8), // CG0-CG3
    ClosureDummy(u8),  // G0-G3

    // For atoms that don't match any of the above
    Unknown,
}

impl AtomType {
    /// Returns the radius of the atom type in Angstroms
    pub fn radius(&self) -> f64 {
        // Radii from AutoDock Vina source (xs_vdw_radii in atom_constants.h)
        match self {
            AtomType::Carbon => 1.9,
            AtomType::Nitrogen => 1.8,
            AtomType::NitrogenH => 1.8,
            AtomType::Oxygen => 1.7,
            AtomType::OxygenH => 1.7,
            AtomType::Sulfur => 2.0,
            AtomType::SulfurH => 2.0,
            AtomType::Phosphorus => 2.1,
            AtomType::Fluorine => 1.5,
            AtomType::Chlorine => 1.8,
            AtomType::Bromine => 2.0,
            AtomType::Iodine => 2.2,
            AtomType::Hydrogen => 1.0,
            AtomType::HydrogenD => 1.0,
            AtomType::Zinc => 1.2,    // Metal donor
            AtomType::Calcium => 1.2, // Metal donor
            AtomType::Manganese => 1.2,
            AtomType::Magnesium => 1.2,
            AtomType::Iron => 1.2,
            // Vina's generic metal radius
            AtomType::Sodium
            | AtomType::Potassium
            | AtomType::Copper
            | AtomType::Nickel
            | AtomType::Cobalt
            | AtomType::Cadmium
            | AtomType::Mercury
            | AtomType::Uranium => 1.2,
            AtomType::ZincPseudo => 0.25,
            AtomType::CarbonClosure(_) => 1.9,
            AtomType::ClosureDummy(_) => 0.0,
            AtomType::Unknown => 1.9, // Default to carbon radius
        }
    }

    /// Covalent radius in Angstroms, which Vina's bond inference uses
    /// (`covalent_radius` in `atom_kind_data`, atom_constants.h). Vina gives
    /// its other metals `metal_covalent_radius`; types Vina does not know get
    /// the same.
    pub fn covalent_radius(&self) -> f64 {
        match self {
            AtomType::Carbon | AtomType::CarbonClosure(_) | AtomType::ClosureDummy(_) => 0.77,
            AtomType::Nitrogen | AtomType::NitrogenH => 0.75,
            AtomType::Oxygen | AtomType::OxygenH => 0.73,
            AtomType::Phosphorus => 1.06,
            AtomType::Sulfur | AtomType::SulfurH => 1.02,
            AtomType::Hydrogen | AtomType::HydrogenD => 0.37,
            AtomType::Fluorine => 0.71,
            AtomType::Iodine => 1.33,
            AtomType::Magnesium => 1.30,
            AtomType::Manganese => 1.39,
            AtomType::Zinc => 1.31,
            AtomType::Calcium => 1.74,
            AtomType::Iron => 1.25,
            AtomType::Chlorine => 0.99,
            AtomType::Bromine => 1.14,
            AtomType::ZincPseudo => 0.0,
            AtomType::Sodium
            | AtomType::Potassium
            | AtomType::Copper
            | AtomType::Nickel
            | AtomType::Cobalt
            | AtomType::Cadmium
            | AtomType::Mercury
            | AtomType::Uranium
            | AtomType::Unknown => 1.75,
        }
    }

    /// Parse atom type from string representation in PDBQT format
    pub fn from_pdbqt_string(s: &str) -> Self {
        // Vina matches type names case-sensitively, which is what tells
        // sodium (`Na`) from an acceptor nitrogen (`NA`). Its non-AD4 metals
        // are matched exactly first; everything else stays case-insensitive.
        match s.trim() {
            "Na" => return AtomType::Sodium,
            "K" => return AtomType::Potassium,
            "Cu" => return AtomType::Copper,
            "Ni" => return AtomType::Nickel,
            "Co" => return AtomType::Cobalt,
            "Cd" => return AtomType::Cadmium,
            "Hg" => return AtomType::Mercury,
            "U" => return AtomType::Uranium,
            _ => {}
        }
        match s.trim().to_uppercase().as_str() {
            "A" => AtomType::Carbon, // Aromatic carbon - same as aliphatic in Vina
            "C" => AtomType::Carbon,
            "N" => AtomType::Nitrogen,
            "NA" => AtomType::NitrogenH, // Nitrogen H-bond acceptor
            "NS" => AtomType::NitrogenH, // Nitrogen in planar system
            "O" => AtomType::Oxygen,
            "OA" => AtomType::OxygenH, // Oxygen H-bond acceptor
            "S" => AtomType::Sulfur,
            "SA" => AtomType::SulfurH, // Sulfur H-bond acceptor
            "P" => AtomType::Phosphorus,
            "F" => AtomType::Fluorine,
            "CL" => AtomType::Chlorine,
            "BR" => AtomType::Bromine,
            "I" => AtomType::Iodine,
            "H" => AtomType::Hydrogen,
            "HD" => AtomType::HydrogenD, // Hydrogen H-bond donor
            "ZN" => AtomType::Zinc,
            "CA" => AtomType::Calcium,
            "MN" => AtomType::Manganese,
            "MG" => AtomType::Magnesium,
            "FE" => AtomType::Iron,
            "TZ" => AtomType::ZincPseudo,
            "CG0" => AtomType::CarbonClosure(0),
            "CG1" => AtomType::CarbonClosure(1),
            "CG2" => AtomType::CarbonClosure(2),
            "CG3" => AtomType::CarbonClosure(3),
            "G0" => AtomType::ClosureDummy(0),
            "G1" => AtomType::ClosureDummy(1),
            "G2" => AtomType::ClosureDummy(2),
            "G3" => AtomType::ClosureDummy(3),
            _ => AtomType::Unknown,
        }
    }

    /// Whether this element can be hydrophobic at all under Vina's `xs` typing.
    ///
    /// Carbon and the halogens qualify. Carbon additionally requires that it not
    /// be bonded to a heteroatom (`C_H` vs `C_P`), which this type-only check
    /// cannot see — see [`Atom::is_hydrophobic`].
    pub fn is_potentially_hydrophobic(&self) -> bool {
        matches!(
            self,
            AtomType::Carbon
                | AtomType::Fluorine
                | AtomType::Chlorine
                | AtomType::Bromine
                | AtomType::Iodine
        )
    }

    /// Whether this element accepts hydrogen bonds.
    ///
    /// Vina treats only the acceptor-flagged PDBQT types `NA` and `OA` as
    /// acceptors. `SA` is not one ("X-Score formulation apparently ignores
    /// SA", `model::assign_types`), and plain `N`/`O`/`S` are non-polar.
    pub fn is_hbond_acceptor(&self) -> bool {
        matches!(self, AtomType::NitrogenH | AtomType::OxygenH)
    }

    /// Whether this is a nitrogen or oxygen of any PDBQT type, the elements
    /// Vina lets donate a hydrogen bond through a bonded polar hydrogen.
    pub fn is_nitrogen_or_oxygen(&self) -> bool {
        matches!(
            self,
            AtomType::Nitrogen | AtomType::NitrogenH | AtomType::Oxygen | AtomType::OxygenH
        )
    }

    /// Whether this is a metal. Vina types every metal `Met_D`: a hydrogen-bond
    /// donor and a heteroatom.
    pub fn is_metal(&self) -> bool {
        matches!(
            self,
            AtomType::Zinc
                | AtomType::Calcium
                | AtomType::Manganese
                | AtomType::Magnesium
                | AtomType::Iron
                | AtomType::Sodium
                | AtomType::Potassium
                | AtomType::Copper
                | AtomType::Nickel
                | AtomType::Cobalt
                | AtomType::Cadmium
                | AtomType::Mercury
                | AtomType::Uranium
        )
    }

    /// Whether this is a hydrogen, which Vina's heavy-atom scoring ignores.
    pub fn is_hydrogen(&self) -> bool {
        matches!(self, AtomType::Hydrogen | AtomType::HydrogenD)
    }

    /// Whether this is a macrocycle closure dummy (`G0`-`G3`).
    pub fn is_closure_dummy(&self) -> bool {
        matches!(self, AtomType::ClosureDummy(_))
    }

    /// Whether `a` and `b` are a closure dummy and the carbon it is glued
    /// onto (`G0` with `CG0`, and so on). Vina pulls such a pair together
    /// during the search to close the ring, and puts it at optimal distance
    /// zero.
    pub fn is_glue_pair(a: AtomType, b: AtomType) -> bool {
        matches!(
            (a, b),
            (AtomType::CarbonClosure(i), AtomType::ClosureDummy(j))
                | (AtomType::ClosureDummy(j), AtomType::CarbonClosure(i)) if i == j
        )
    }

    /// Whether Vina's pair potentials score `a` against `b` at all: never a
    /// hydrogen, and a closure dummy only against its own closure carbon
    /// (Vina's `is_unmatched_closure_dummy`, and its grids skip dummies).
    pub fn pair_is_scored(a: AtomType, b: AtomType) -> bool {
        if a.is_hydrogen() || b.is_hydrogen() {
            return false;
        }
        if a.is_closure_dummy() || b.is_closure_dummy() {
            return Self::is_glue_pair(a, b);
        }
        true
    }

    /// Whether this is a heteroatom for the purpose of the `C_H`/`C_P` split.
    ///
    /// As Vina's `ad_is_heteroatom` plus metals: every type but carbon and
    /// hydrogen, so a carbon bonded to a halogen or a metal is polar too.
    /// That includes the closure types `CG`/`G`, which Vina does not count
    /// as carbon either.
    pub fn is_heteroatom(&self) -> bool {
        !matches!(
            self,
            AtomType::Carbon
                | AtomType::Hydrogen
                | AtomType::HydrogenD
                | AtomType::ZincPseudo
                | AtomType::Unknown
        )
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
            AtomType::Sodium => "Na",
            AtomType::Potassium => "K",
            AtomType::Copper => "Cu",
            AtomType::Nickel => "Ni",
            AtomType::Cobalt => "Co",
            AtomType::Cadmium => "Cd",
            AtomType::Mercury => "Hg",
            AtomType::Uranium => "U",
            AtomType::ZincPseudo => "TZ",
            AtomType::CarbonClosure(n) => ["CG0", "CG1", "CG2", "CG3"][(n & 3) as usize],
            AtomType::ClosureDummy(n) => ["G0", "G1", "G2", "G3"][(n & 3) as usize],
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
    pub residue_num: i32,

    /// Chain identifier
    pub chain_id: char,

    /// Partial charge
    pub charge: f64,

    /// Is this atom part of a flexible side chain?
    pub is_flexible: bool,

    /// Whether this atom is hydrophobic under Vina's `xs` typing.
    ///
    /// Carbon and the halogens qualify, but carbon only when it is *not* bonded
    /// to a heteroatom — Vina's `C_H` vs `C_P` distinction. That depends on the
    /// bond graph rather than the atom type alone, so this is assigned by
    /// [`Molecule::assign_xs_types`](crate::molecule::Molecule::assign_xs_types)
    /// once connectivity is known. The value here is the type-only
    /// approximation until then.
    pub is_hydrophobic: bool,

    /// Whether this atom donates a hydrogen bond, i.e. it is an N/O/S carrying
    /// a bonded polar hydrogen.
    ///
    /// Also bond-graph dependent, and assigned by the same pass. A receptor
    /// prepared without polar hydrogens therefore has no donors at all, which
    /// is what Vina does too.
    pub is_hbond_donor: bool,
}

impl Atom {
    /// Create a new atom
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        atom_type: AtomType,
        coordinates: Vector3<f64>,
        name: String,
        serial: u32,
        residue_name: String,
        residue_num: i32,
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
            is_hydrophobic: atom_type.is_potentially_hydrophobic(),
            is_hbond_donor: false,
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

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_atom_type_radius() {
        // Radii from AutoDock Vina source
        assert_eq!(AtomType::Carbon.radius(), 1.9);
        assert_eq!(AtomType::Nitrogen.radius(), 1.8);
        assert_eq!(AtomType::Oxygen.radius(), 1.7);
        assert_eq!(AtomType::Hydrogen.radius(), 1.0);
        assert_eq!(AtomType::ZincPseudo.radius(), 0.25);
    }

    #[test]
    fn test_atom_type_from_pdbqt_string() {
        assert_eq!(AtomType::from_pdbqt_string("C"), AtomType::Carbon);
        assert_eq!(AtomType::from_pdbqt_string("NA"), AtomType::NitrogenH);
        assert_eq!(AtomType::from_pdbqt_string("OA"), AtomType::OxygenH);
        assert_eq!(AtomType::from_pdbqt_string("HD"), AtomType::HydrogenD);
        assert_eq!(AtomType::from_pdbqt_string("TZ"), AtomType::ZincPseudo);
        assert_eq!(AtomType::from_pdbqt_string("UNKNOWN"), AtomType::Unknown);
    }

    #[test]
    fn test_closure_types() {
        assert_eq!(
            AtomType::from_pdbqt_string("CG2"),
            AtomType::CarbonClosure(2)
        );
        assert_eq!(AtomType::from_pdbqt_string("G3"), AtomType::ClosureDummy(3));
        assert_eq!(AtomType::CarbonClosure(1).to_pdbqt_string(), "CG1");
        assert_eq!(AtomType::ClosureDummy(0).to_pdbqt_string(), "G0");
        assert!(AtomType::CarbonClosure(0).is_heteroatom());
        assert!(!AtomType::CarbonClosure(0).is_potentially_hydrophobic());

        let (cg0, g0, g1) = (
            AtomType::CarbonClosure(0),
            AtomType::ClosureDummy(0),
            AtomType::ClosureDummy(1),
        );
        assert!(AtomType::is_glue_pair(g0, cg0) && AtomType::pair_is_scored(cg0, g0));
        assert!(!AtomType::is_glue_pair(g1, cg0) && !AtomType::pair_is_scored(g1, cg0));
        assert!(!AtomType::pair_is_scored(g0, AtomType::Carbon));
        assert!(AtomType::pair_is_scored(cg0, AtomType::OxygenH));
    }

    #[test]
    fn test_metal_names_are_case_sensitive() {
        // Open Babel writes sodium as `Na`; upper-casing it made it an
        // acceptor nitrogen. Vina matches type names exactly.
        assert_eq!(AtomType::from_pdbqt_string("Na"), AtomType::Sodium);
        assert_eq!(AtomType::from_pdbqt_string("NA"), AtomType::NitrogenH);
        assert_eq!(AtomType::from_pdbqt_string("K"), AtomType::Potassium);
        assert_eq!(AtomType::from_pdbqt_string("Cu"), AtomType::Copper);
        assert_eq!(AtomType::from_pdbqt_string("Ca"), AtomType::Calcium);
        for metal in [AtomType::Sodium, AtomType::Copper, AtomType::Zinc] {
            assert!(metal.is_metal());
            assert_eq!(AtomType::from_pdbqt_string(metal.to_pdbqt_string()), metal);
        }
    }

    #[test]
    fn test_vina_acceptors_are_na_and_oa_only() {
        assert!(AtomType::NitrogenH.is_hbond_acceptor());
        assert!(AtomType::OxygenH.is_hbond_acceptor());
        assert!(!AtomType::SulfurH.is_hbond_acceptor());
        assert!(!AtomType::Nitrogen.is_hbond_acceptor());
    }

    #[test]
    fn test_atom_type_to_pdbqt_string() {
        assert_eq!(AtomType::Carbon.to_pdbqt_string(), "C");
        assert_eq!(AtomType::NitrogenH.to_pdbqt_string(), "NA");
        assert_eq!(AtomType::OxygenH.to_pdbqt_string(), "OA");
        assert_eq!(AtomType::HydrogenD.to_pdbqt_string(), "HD");
        assert_eq!(AtomType::ZincPseudo.to_pdbqt_string(), "TZ");
        assert_eq!(AtomType::Unknown.to_pdbqt_string(), "X");
    }

    #[test]
    fn test_atom_creation() {
        let atom = Atom::new(
            AtomType::Carbon,
            Vector3::new(1.0, 2.0, 3.0),
            "CA".to_string(),
            1,
            "ALA".to_string(),
            1,
            'A',
            0.0,
        );

        assert_eq!(atom.atom_type, AtomType::Carbon);
        assert_eq!(atom.coordinates, Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(atom.name, "CA");
        assert_eq!(atom.serial, 1);
        assert_eq!(atom.residue_name, "ALA");
        assert_eq!(atom.residue_num, 1);
        assert_eq!(atom.chain_id, 'A');
        assert_eq!(atom.charge, 0.0);
        assert!(!atom.is_flexible);
    }

    #[test]
    fn test_atom_distance() {
        let atom1 = Atom::new(
            AtomType::Carbon,
            Vector3::new(0.0, 0.0, 0.0),
            "CA".to_string(),
            1,
            "ALA".to_string(),
            1,
            'A',
            0.0,
        );

        let atom2 = Atom::new(
            AtomType::Carbon,
            Vector3::new(1.0, 1.0, 1.0),
            "CA".to_string(),
            2,
            "ALA".to_string(),
            1,
            'A',
            0.0,
        );

        // Distance should be sqrt(3) ≈ 1.732
        assert!((atom1.distance(&atom2) - 1.732).abs() < 0.001);
    }

    #[test]
    fn test_h_bond_donor() {
        let h_donor = Atom::new(
            AtomType::HydrogenD,
            Vector3::new(0.0, 0.0, 0.0),
            "H".to_string(),
            1,
            "ALA".to_string(),
            1,
            'A',
            0.0,
        );

        let non_donor = Atom::new(
            AtomType::Carbon,
            Vector3::new(0.0, 0.0, 0.0),
            "C".to_string(),
            2,
            "ALA".to_string(),
            1,
            'A',
            0.0,
        );

        assert!(h_donor.is_h_bond_donor());
        assert!(!non_donor.is_h_bond_donor());
    }

    #[test]
    fn test_h_bond_acceptor() {
        let o_acceptor = Atom::new(
            AtomType::OxygenH,
            Vector3::new(0.0, 0.0, 0.0),
            "O".to_string(),
            1,
            "ALA".to_string(),
            1,
            'A',
            0.0,
        );

        let n_acceptor = Atom::new(
            AtomType::NitrogenH,
            Vector3::new(0.0, 0.0, 0.0),
            "N".to_string(),
            2,
            "ALA".to_string(),
            1,
            'A',
            0.0,
        );

        let non_acceptor = Atom::new(
            AtomType::Carbon,
            Vector3::new(0.0, 0.0, 0.0),
            "C".to_string(),
            3,
            "ALA".to_string(),
            1,
            'A',
            0.0,
        );

        assert!(o_acceptor.is_h_bond_acceptor());
        assert!(n_acceptor.is_h_bond_acceptor());
        assert!(!non_acceptor.is_h_bond_acceptor());
    }

    #[test]
    fn test_atom_display() {
        let atom = Atom::new(
            AtomType::Carbon,
            Vector3::new(1.0, 2.0, 3.0),
            "CA".to_string(),
            1,
            "ALA".to_string(),
            1,
            'A',
            0.5,
        );

        let expected = "C(1, 2, 3) [0.5]";
        assert_eq!(format!("{}", atom), expected);
    }
}
