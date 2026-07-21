//! High-level scoring of a ligand pose with the Vina force field.

use crate::forcefield::vina::VinaForceField;
use crate::forcefield::ForceField;
use crate::molecule::Molecule;

/// Scores a ligand pose against a receptor using the Vina scoring function.
///
/// This wraps the same energy terms the optimizer uses, so a pose scored here
/// and the same pose produced by docking agree.
///
/// ```no_run
/// use rustdock_vina::{parse_pdbqt, VinaForceField, VinaScore};
///
/// # fn main() -> rustdock_vina::Result<()> {
/// let receptor = parse_pdbqt("receptor.pdbqt")?;
/// let ligand = parse_pdbqt("ligand.pdbqt")?;
/// let forcefield = VinaForceField::default();
///
/// let score = VinaScore::new(&ligand, &forcefield);
/// println!("affinity: {:.3} kcal/mol", score.affinity(&receptor)?);
/// # Ok(())
/// # }
/// ```
pub struct VinaScore<'a> {
    molecule: &'a Molecule,
    forcefield: &'a VinaForceField,
}

impl<'a> VinaScore<'a> {
    /// Create a new Vina scoring function
    pub fn new(molecule: &'a Molecule, forcefield: &'a VinaForceField) -> Self {
        Self {
            molecule,
            forcefield,
        }
    }

    /// Intermolecular interaction energy with `receptor`, before the torsional
    /// normalization.
    pub fn intermolecular_energy(&self, receptor: &Molecule) -> Result<f64, crate::Error> {
        let mut total = 0.0;
        for lig_atom in &self.molecule.atoms {
            for rec_atom in &receptor.atoms {
                let d = lig_atom.distance(rec_atom);
                if d > 0.01 {
                    total += self.forcefield.atom_pair_energy(lig_atom, rec_atom, d)?;
                }
            }
        }
        Ok(total)
    }

    /// The ligand's internal energy, over pairs separated by more than three
    /// bonds.
    ///
    /// Vina reports this alongside the affinity but does not include it there:
    /// it cancels against the unbound reference state.
    pub fn internal_energy(&self) -> Result<f64, crate::Error> {
        let mut total = 0.0;
        for (i, j) in self.molecule.intramolecular_pairs() {
            let (a, b) = (&self.molecule.atoms[i], &self.molecule.atoms[j]);
            let d = a.distance(b);
            if d > 0.01 {
                total += self.forcefield.atom_pair_energy(a, b, d)?;
            }
        }
        Ok(total)
    }

    /// Predicted binding affinity in kcal/mol.
    ///
    /// This is the intermolecular energy divided by `1 + w_rot * N_rot`, the
    /// same quantity Vina prints as `VINA RESULT`.
    pub fn affinity(&self, receptor: &Molecule) -> Result<f64, crate::Error> {
        let inter = self.intermolecular_energy(receptor)?;
        Ok(self
            .forcefield
            .calculate_affinity(inter, self.molecule.torsions.len()))
    }
}
