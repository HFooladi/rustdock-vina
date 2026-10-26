pub mod vina_score;

pub use vina_score::VinaScore;

use crate::forcefield::ForceField;
use crate::molecule::Molecule;

/// Pairs closer than this are skipped as overlapping atoms.
const MIN_DISTANCE: f64 = 0.01;

/// Pairs farther than this contribute nothing with Vina or Vinardo.
const CUTOFF: f64 = 8.0;

/// The score of a ligand pose against a receptor, as `rustdock score` prints it.
#[derive(Debug, Clone, PartialEq)]
pub struct PoseScore {
    /// Sum of the force field's pair energies, before torsional normalization.
    pub intermolecular: f64,
    /// Predicted binding affinity in kcal/mol: `intermolecular` after the
    /// force field's torsional normalization.
    pub affinity: f64,
    /// Vina's weighted torsion count of the ligand.
    pub num_tors: f64,
    /// Number of ligand–receptor pairs within the cutoff.
    pub pairs: usize,
    /// Sum of the repulsive (positive) pair energies.
    pub positive: f64,
    /// Sum of the attractive (negative) pair energies.
    pub negative: f64,
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
    let mut intermolecular = 0.0;
    let mut pairs = 0;
    let mut positive = 0.0;
    let mut negative = 0.0;
    for lig_atom in &ligand.atoms {
        for rec_atom in &receptor.atoms {
            let distance = lig_atom.distance(rec_atom);
            if distance < CUTOFF && distance > MIN_DISTANCE {
                let energy = forcefield.atom_pair_energy(lig_atom, rec_atom, distance)?;
                intermolecular += energy;
                pairs += 1;
                if energy > 0.0 {
                    positive += energy;
                } else {
                    negative += energy;
                }
            }
        }
    }
    let num_tors = ligand.vina_torsion_count();
    let affinity = forcefield.apply_torsional_normalization(intermolecular, num_tors);
    Ok(PoseScore {
        intermolecular,
        affinity,
        num_tors,
        pairs,
        positive,
        negative,
    })
}
