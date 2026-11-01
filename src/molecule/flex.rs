//! Flexible receptor side chains, as Vina reads them with `--flex`.

use super::Molecule;
use crate::atom::Atom;
use std::collections::VecDeque;

/// The flexible side chains of a receptor: every `BEGIN_RES`…`END_RES` block
/// of a flex PDBQT, read by [`parse_flex_pdbqt`](crate::io::parse_flex_pdbqt).
///
/// Each BRANCH is a torsion. A residue's ROOT atoms (Meeko and ADT write the
/// CA) stay fixed, and so does the first atom of each branch off the ROOT
/// (the CB), which lies on the branch's rotation axis. As in Vina, the fixed atoms are part of neither
/// the rigid receptor nor the ligand's partners: they only interact with
/// the residues' moving atoms. The rigid receptor file must leave all of
/// these atoms out, as Meeko and ADT write it.
#[derive(Debug, Clone)]
pub struct FlexResidues {
    /// The atoms of all residues in file order, with their bonds and
    /// torsions (one per BRANCH), and the file's records for writing.
    pub molecule: Molecule,
    /// Whether each atom is fixed: a ROOT atom, or the first atom of a
    /// branch off the ROOT.
    pub fixed: Vec<bool>,
    /// Index of each atom's residue, in file order.
    pub residue_of: Vec<usize>,
    /// Rigid body of each atom, as Vina's torsion tree has them: the ROOT
    /// is body 0 and each BRANCH is a body, except that a branch's first
    /// atom, on its rotation axis, belongs to the parent's body.
    pub(crate) body_of: Vec<usize>,
    /// The two atoms of each body's rotation axis: the atom it hangs from
    /// and the branch's first atom (none for the ROOT).
    pub(crate) axis: Vec<(Option<usize>, Option<usize>)>,
}

impl FlexResidues {
    /// Number of residues.
    pub fn num_residues(&self) -> usize {
        self.residue_of.last().map_or(0, |&r| r + 1)
    }

    /// Indices of the atoms that move with the torsions.
    pub fn movable(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.fixed.len()).filter(|&i| !self.fixed[i])
    }

    /// The moving atoms, which score against the rigid receptor and the
    /// ligand.
    pub fn movable_atoms(&self) -> impl Iterator<Item = &Atom> + '_ {
        self.movable().map(|i| &self.molecule.atoms[i])
    }

    /// Whether atoms `a` and `b` can move relative to each other (Vina's
    /// `DISTANCE_VARIABLE`): not both fixed, in different rigid bodies, and
    /// neither on the axis of the other's body.
    pub(crate) fn variable(&self, a: usize, b: usize) -> bool {
        let (ba, bb) = (self.body_of[a], self.body_of[b]);
        let on_axis = |atom: usize, body: usize| {
            let (begin, end) = self.axis[body];
            begin == Some(atom) || end == Some(atom)
        };
        !(self.fixed[a] && self.fixed[b]) && ba != bb && !on_axis(a, bb) && !on_axis(b, ba)
    }

    /// The pairs Vina scores within and between the residues (its
    /// `other_pairs`): heavy atoms that can move relative to each other and
    /// are more than three bonds apart. Fixed atoms take part, but two fixed
    /// atoms never move relative to each other.
    pub fn internal_pairs(&self) -> Vec<(usize, usize)> {
        let atoms = &self.molecule.atoms;
        let n = atoms.len();
        let mut adjacency = vec![Vec::new(); n];
        for bond in &self.molecule.bonds {
            if bond.atom1_idx < n && bond.atom2_idx < n {
                adjacency[bond.atom1_idx].push(bond.atom2_idx);
                adjacency[bond.atom2_idx].push(bond.atom1_idx);
            }
        }

        let mut pairs = Vec::new();
        let mut depth = vec![usize::MAX; n];
        let mut queue = VecDeque::new();
        for i in 0..n {
            if atoms[i].atom_type.is_hydrogen() {
                continue;
            }
            depth.iter_mut().for_each(|d| *d = usize::MAX);
            depth[i] = 0;
            queue.clear();
            queue.push_back(i);
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
            for j in (i + 1)..n {
                if depth[j] > 3 && !atoms[j].atom_type.is_hydrogen() && self.variable(i, j) {
                    pairs.push((i, j));
                }
            }
        }
        pairs
    }
}
