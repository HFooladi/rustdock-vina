//! Molecule representation and related functionality

use crate::atom::{Atom, AtomType};
use nalgebra::{Unit, UnitQuaternion, Vector3};
use std::collections::HashMap;
use std::f64::consts::PI;
use thiserror::Error;

/// Wrap an angle into `[-PI, PI]`.
fn wrap_angle(mut angle: f64) -> f64 {
    while angle > PI {
        angle -= 2.0 * PI;
    }
    while angle < -PI {
        angle += 2.0 * PI;
    }
    angle
}

/// Errors that can occur when working with molecules
#[derive(Error, Debug)]
pub enum MoleculeError {
    #[error("Invalid bond: atoms {0} and {1} not found")]
    InvalidBond(u32, u32),

    #[error("Invalid atom index: {0}")]
    InvalidAtomIndex(usize),

    #[error("No atoms in molecule")]
    EmptyMolecule,
}

/// Represents a chemical bond between two atoms
#[derive(Debug, Clone)]
pub struct Bond {
    /// Index of the first atom
    pub atom1_idx: usize,

    /// Index of the second atom
    pub atom2_idx: usize,

    /// Is this bond rotatable?
    pub rotatable: bool,
}

/// Represents a molecular torsion (rotatable bond and the atoms affected by rotation)
#[derive(Debug, Clone)]
pub struct Torsion {
    /// Index of the rotatable bond
    pub bond_idx: usize,

    /// Indices of atoms that move when this bond rotates
    pub moving_atoms: Vec<usize>,

    /// Current torsion angle (in radians)
    pub angle: f64,
}

/// Represents a molecule (ligand or receptor)
#[derive(Debug, Clone)]
pub struct Molecule {
    /// Name of the molecule
    pub name: String,

    /// List of atoms in the molecule
    pub atoms: Vec<Atom>,

    /// List of bonds between atoms
    pub bonds: Vec<Bond>,

    /// List of rotatable bonds (torsions)
    pub torsions: Vec<Torsion>,

    /// Is this molecule rigid or has flexible parts?
    pub is_rigid: bool,
}

impl Molecule {
    /// Create a new empty molecule
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            atoms: Vec::new(),
            bonds: Vec::new(),
            torsions: Vec::new(),
            is_rigid: true,
        }
    }

    /// Add an atom to the molecule
    pub fn add_atom(&mut self, atom: Atom) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(atom);
        idx
    }

    /// Add a bond between two atoms
    pub fn add_bond(
        &mut self,
        atom1_idx: usize,
        atom2_idx: usize,
        rotatable: bool,
    ) -> Result<usize, MoleculeError> {
        if atom1_idx >= self.atoms.len() || atom2_idx >= self.atoms.len() {
            return Err(MoleculeError::InvalidBond(
                atom1_idx as u32,
                atom2_idx as u32,
            ));
        }

        let bond = Bond {
            atom1_idx,
            atom2_idx,
            rotatable,
        };

        let idx = self.bonds.len();
        self.bonds.push(bond);

        if rotatable {
            self.is_rigid = false;
        }

        Ok(idx)
    }

    /// Assign the bond-graph-dependent parts of Vina's `xs` atom typing.
    ///
    /// Two properties cannot be derived from an atom's element alone, and
    /// getting them wrong systematically inflates the score:
    ///
    /// - **Hydrophobicity.** Vina splits carbon into `C_H` (bonded only to
    ///   carbon and hydrogen) and `C_P` (bonded to a heteroatom); only `C_H`
    ///   is hydrophobic. Treating every carbon as hydrophobic over-rewards
    ///   buried polar groups — in a protein roughly half of all carbons are
    ///   `C_P`, since each residue contributes a backbone CA and carbonyl C.
    /// - **Donor status.** An N/O/S donates a hydrogen bond only when it
    ///   actually carries a polar hydrogen. Inferring it from the acceptor
    ///   labels `NA`/`OA`/`SA` instead fabricates donors, and in a structure
    ///   prepared without polar hydrogens it fabricates every one of them.
    ///
    /// Call this after connectivity is known; [`parse_pdbqt`](crate::io::parse_pdbqt)
    /// does so automatically.
    pub fn assign_xs_types(&mut self) {
        let n = self.atoms.len();
        let mut bonded_to_hetero = vec![false; n];
        let mut bonded_to_polar_h = vec![false; n];

        for bond in &self.bonds {
            let (i, j) = (bond.atom1_idx, bond.atom2_idx);
            if i >= n || j >= n {
                continue;
            }
            let (ti, tj) = (self.atoms[i].atom_type, self.atoms[j].atom_type);

            if tj.is_heteroatom() {
                bonded_to_hetero[i] = true;
            }
            if ti.is_heteroatom() {
                bonded_to_hetero[j] = true;
            }
            // Only HD, the PDBQT polar hydrogen, makes its partner a donor.
            if tj == AtomType::HydrogenD {
                bonded_to_polar_h[i] = true;
            }
            if ti == AtomType::HydrogenD {
                bonded_to_polar_h[j] = true;
            }
        }

        for (i, atom) in self.atoms.iter_mut().enumerate() {
            let t = atom.atom_type;
            atom.is_hydrophobic =
                t.is_potentially_hydrophobic() && !(t == AtomType::Carbon && bonded_to_hetero[i]);
            atom.is_hbond_donor = t.is_hbond_acceptor() && bonded_to_polar_h[i];
        }
    }

    /// Atom index pairs that contribute to the intramolecular energy.
    ///
    /// Vina scores a ligand against itself over pairs separated by more than
    /// three bonds; closer pairs are held rigid by the covalent geometry and
    /// would contribute a large constant repulsion. Hydrogens are excluded, as
    /// everywhere else in the scoring.
    ///
    /// The result depends only on topology, not coordinates, so it is safe to
    /// compute once and reuse across an entire optimization.
    pub fn intramolecular_pairs(&self) -> Vec<(usize, usize)> {
        let n = self.atoms.len();
        let mut adjacency = vec![Vec::new(); n];
        for bond in &self.bonds {
            if bond.atom1_idx < n && bond.atom2_idx < n {
                adjacency[bond.atom1_idx].push(bond.atom2_idx);
                adjacency[bond.atom2_idx].push(bond.atom1_idx);
            }
        }

        let mut pairs = Vec::new();
        let mut depth = vec![usize::MAX; n];
        let mut queue = std::collections::VecDeque::new();

        for start in 0..n {
            if self.atoms[start].atom_type.is_hydrogen() {
                continue;
            }

            // Breadth-first walk out to three bonds; anything further away (or
            // unreachable, in a disconnected structure) is a scoring pair.
            depth.iter_mut().for_each(|d| *d = usize::MAX);
            depth[start] = 0;
            queue.clear();
            queue.push_back(start);
            while let Some(cur) = queue.pop_front() {
                if depth[cur] >= 3 {
                    continue;
                }
                for &next in &adjacency[cur] {
                    if depth[next] == usize::MAX {
                        depth[next] = depth[cur] + 1;
                        queue.push_back(next);
                    }
                }
            }

            for (other, &bond_distance) in depth.iter().enumerate().skip(start + 1) {
                if self.atoms[other].atom_type.is_hydrogen() {
                    continue;
                }
                if bond_distance > 3 {
                    pairs.push((start, other));
                }
            }
        }

        pairs
    }

    /// Rotate a torsion by `delta_angle` radians, moving the atoms distal to
    /// the rotatable bond.
    ///
    /// This is the single implementation of torsional motion: both the Monte
    /// Carlo moves and the local optimizer's torsion gradients go through it,
    /// so a conformation's torsion angles always match its coordinates.
    pub fn rotate_torsion(
        &mut self,
        torsion_idx: usize,
        delta_angle: f64,
    ) -> Result<(), MoleculeError> {
        let torsion = self
            .torsions
            .get(torsion_idx)
            .ok_or(MoleculeError::InvalidAtomIndex(torsion_idx))?;
        let bond = self
            .bonds
            .get(torsion.bond_idx)
            .ok_or(MoleculeError::InvalidAtomIndex(torsion.bond_idx))?;
        let (a1, a2) = (bond.atom1_idx, bond.atom2_idx);

        if a1 >= self.atoms.len() || a2 >= self.atoms.len() {
            return Err(MoleculeError::InvalidBond(a1 as u32, a2 as u32));
        }
        if torsion.moving_atoms.iter().any(|&i| i >= self.atoms.len()) {
            return Err(MoleculeError::InvalidAtomIndex(self.atoms.len()));
        }

        let origin = self.atoms[a1].coordinates;
        let axis_vec = self.atoms[a2].coordinates - origin;
        if axis_vec.norm() < 1e-10 {
            return Err(MoleculeError::InvalidBond(a1 as u32, a2 as u32));
        }
        let rotation = UnitQuaternion::from_axis_angle(&Unit::new_normalize(axis_vec), delta_angle);

        let moving = self.torsions[torsion_idx].moving_atoms.clone();
        for idx in moving {
            let pos = self.atoms[idx].coordinates - origin;
            self.atoms[idx].coordinates = rotation.transform_vector(&pos) + origin;
        }

        let torsion = &mut self.torsions[torsion_idx];
        torsion.angle = wrap_angle(torsion.angle + delta_angle);
        Ok(())
    }

    /// Translate every atom by `offset`.
    pub fn translate(&mut self, offset: &Vector3<f64>) {
        for atom in &mut self.atoms {
            atom.coordinates += offset;
        }
    }

    /// Rotate the whole molecule about its centroid.
    pub fn rotate_about_center(&mut self, rotation: &UnitQuaternion<f64>) {
        let Ok(center) = self.center() else { return };
        for atom in &mut self.atoms {
            let pos = atom.coordinates - center;
            atom.coordinates = rotation.transform_vector(&pos) + center;
        }
    }

    /// Whether every atom lies inside the axis-aligned search box.
    ///
    /// Note this is a *stricter* test than docking should use — see
    /// [`center_within_box`](Self::center_within_box). A drug-like ligand is
    /// often longer than the box: imatinib spans 21 A, so it never fits
    /// entirely inside the conventional 20 A box at any orientation.
    pub fn is_within_box(&self, center: &Vector3<f64>, box_size: &Vector3<f64>) -> bool {
        let half = box_size * 0.5;
        let (lo, hi) = (center - half, center + half);
        self.atoms.iter().all(|atom| {
            let p = &atom.coordinates;
            p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y && p.z >= lo.z && p.z <= hi.z
        })
    }

    /// Whether the molecule's centroid lies inside the search box.
    ///
    /// This is the containment rule the search uses, and it is what Vina
    /// enforces: the box bounds where the ligand may be *placed*, and atoms are
    /// allowed to extend beyond it. Requiring total containment instead
    /// rejects every pose of any ligand larger than the box, which silently
    /// turns the search into a no-op.
    pub fn center_within_box(&self, center: &Vector3<f64>, box_size: &Vector3<f64>) -> bool {
        let Ok(c) = self.center() else { return false };
        let half = box_size * 0.5;
        (c.x - center.x).abs() <= half.x
            && (c.y - center.y).abs() <= half.y
            && (c.z - center.z).abs() <= half.z
    }

    /// Get the center of the molecule
    pub fn center(&self) -> Result<Vector3<f64>, MoleculeError> {
        if self.atoms.is_empty() {
            return Err(MoleculeError::EmptyMolecule);
        }

        let sum = self
            .atoms
            .iter()
            .fold(Vector3::zeros(), |acc, atom| acc + atom.coordinates);

        Ok(sum / self.atoms.len() as f64)
    }

    /// Get the bounding box of the molecule
    pub fn bounding_box(&self) -> Result<(Vector3<f64>, Vector3<f64>), MoleculeError> {
        if self.atoms.is_empty() {
            return Err(MoleculeError::EmptyMolecule);
        }

        let mut min = Vector3::new(f64::MAX, f64::MAX, f64::MAX);
        let mut max = Vector3::new(f64::MIN, f64::MIN, f64::MIN);

        for atom in &self.atoms {
            // Update minimum corner
            min.x = min.x.min(atom.coordinates.x);
            min.y = min.y.min(atom.coordinates.y);
            min.z = min.z.min(atom.coordinates.z);

            // Update maximum corner
            max.x = max.x.max(atom.coordinates.x);
            max.y = max.y.max(atom.coordinates.y);
            max.z = max.z.max(atom.coordinates.z);
        }

        Ok((min, max))
    }

    /// Count atoms by type
    pub fn count_atom_types(&self) -> HashMap<AtomType, usize> {
        let mut counts = HashMap::new();

        for atom in &self.atoms {
            *counts.entry(atom.atom_type).or_insert(0) += 1;
        }

        counts
    }
}
