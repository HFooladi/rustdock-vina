pub mod vina_score;

pub use vina_score::VinaScore;

use crate::atom::Atom;
use crate::forcefield::ForceField;
use crate::molecule::{FlexResidues, Molecule};

/// Pairs closer than this are skipped as overlapping atoms.
const MIN_DISTANCE: f64 = 0.01;

/// Pairs farther than this contribute nothing with Vina or Vinardo.
const CUTOFF: f64 = 8.0;

/// The score of a ligand pose against a receptor, as `rustdock score` prints it.
#[derive(Debug, Clone, PartialEq)]
pub struct PoseScore {
    /// Sum of the force field's pair energies between the ligand and the
    /// receptor, flexible side chains included, before torsional
    /// normalization.
    pub intermolecular: f64,
    /// Predicted binding affinity in kcal/mol: `intermolecular` after the
    /// force field's torsional normalization.
    pub affinity: f64,
    /// Vina's weighted torsion count of the ligand.
    pub num_tors: f64,
    /// Number of ligand–receptor pairs within the cutoff, flexible side
    /// chains included.
    pub pairs: usize,
    /// Sum of the repulsive (positive) pair energies.
    pub positive: f64,
    /// Sum of the attractive (negative) pair energies.
    pub negative: f64,
    /// The flexible side chains' part of `intermolecular`; zero without them.
    pub ligand_flex: f64,
    /// Flexible side chains against the rigid receptor. Like the next term,
    /// it is part of Vina's intramolecular energy, not of the affinity.
    pub flex_receptor: f64,
    /// Flexible side chains against themselves and each other.
    pub flex_flex: f64,
}

/// Score a ligand pose against a receptor with any force field.
///
/// This is the quantity Vina's `--score_only` reports as the intermolecular
/// energy and the estimated free energy of binding.
pub fn score_pose(
    ligand: &Molecule,
    receptor: &Molecule,
    forcefield: &dyn ForceField,
) -> Result<PoseScore, crate::Error> {
    score_pose_with_flex(ligand, receptor, None, forcefield)
}

/// Score a ligand pose against a receptor with flexible side chains, as
/// `vina --score_only --flex` does; `receptor` is the rigid part.
///
/// The affinity is the normalized intermolecular energy. Vina's own
/// `--score_only` affinity differs by its flex–receptor energy on the grid
/// minus the exact one, which do not cancel; with `--unbound_energy` set to
/// its printed internal energy, Vina reports this affinity.
pub fn score_pose_with_flex(
    ligand: &Molecule,
    receptor: &Molecule,
    flex: Option<&FlexResidues>,
    forcefield: &dyn ForceField,
) -> Result<PoseScore, crate::Error> {
    let mut sum = PairSum::default();
    for lig_atom in &ligand.atoms {
        for rec_atom in &receptor.atoms {
            sum.add(forcefield, lig_atom, rec_atom)?;
        }
    }
    let rigid = sum.energy;

    let (mut flex_receptor, mut flex_flex) = (0.0, 0.0);
    if let Some(flex) = flex {
        for lig_atom in &ligand.atoms {
            for flex_atom in flex.movable_atoms() {
                sum.add(forcefield, lig_atom, flex_atom)?;
            }
        }
        let mut other = PairSum::default();
        for flex_atom in flex.movable_atoms() {
            for rec_atom in &receptor.atoms {
                other.add(forcefield, flex_atom, rec_atom)?;
            }
        }
        flex_receptor = other.energy;
        let atoms = &flex.molecule.atoms;
        let mut other = PairSum::default();
        for (i, j) in flex.internal_pairs() {
            other.add(forcefield, &atoms[i], &atoms[j])?;
        }
        flex_flex = other.energy;
    }

    let intermolecular = sum.energy;
    let num_tors = ligand.vina_torsion_count();
    let affinity = forcefield.apply_torsional_normalization(intermolecular, num_tors);
    Ok(PoseScore {
        intermolecular,
        affinity,
        num_tors,
        pairs: sum.pairs,
        positive: sum.positive,
        negative: sum.negative,
        ligand_flex: intermolecular - rigid,
        flex_receptor,
        flex_flex,
    })
}

/// A running sum of pair energies within the cutoff.
#[derive(Default)]
struct PairSum {
    energy: f64,
    pairs: usize,
    positive: f64,
    negative: f64,
}

impl PairSum {
    fn add(&mut self, forcefield: &dyn ForceField, a: &Atom, b: &Atom) -> Result<(), crate::Error> {
        let distance = a.distance(b);
        if distance < CUTOFF && distance > MIN_DISTANCE {
            let energy = forcefield.atom_pair_energy(a, b, distance)?;
            self.energy += energy;
            self.pairs += 1;
            if energy > 0.0 {
                self.positive += energy;
            } else {
                self.negative += energy;
            }
        }
        Ok(())
    }
}
