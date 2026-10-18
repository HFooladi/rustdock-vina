//! Input/output functionality for molecular docking

use nalgebra::Vector3;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use thiserror::Error;

use crate::atom::{Atom, AtomType};
use crate::molecule::{Bond, Molecule, Torsion};

/// Errors that can occur during file I/O operations
#[derive(Error, Debug)]
pub enum IoError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Parse error at line {line}: {message}")]
    Parse { line: usize, message: String },

    #[error("Invalid file format: {0}")]
    InvalidFormat(String),

    #[error("Unsupported file format: {0}")]
    UnsupportedFormat(String),
}

/// Parse a PDBQT file into a Molecule
pub fn parse_pdbqt<P: AsRef<Path>>(path: P) -> Result<Molecule, IoError> {
    let file = File::open(path.as_ref())?;
    let reader = BufReader::new(file);

    let mut molecule = Molecule::new(
        path.as_ref()
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown"),
    );

    // Keep track of branch structure for building torsions
    let mut branch_stack: Vec<(usize, usize)> = Vec::new();
    let mut branch_atoms: Vec<Vec<usize>> = Vec::new();
    let mut current_branch_atoms: Vec<usize> = Vec::new();
    // Rigid fragment of each atom: the ROOT is 0, each BRANCH opens a new one.
    let mut fragment_of: Vec<usize> = Vec::new();
    let mut fragment_stack: Vec<usize> = vec![0];
    let mut fragment_count = 1;

    // Track model number for multi-model files
    let mut model_number = 1;
    let mut line_number = 0;

    for line in reader.lines() {
        let line = line?;
        line_number += 1;

        // Skip empty lines and comments
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }

        // Keep the file's own records for the writer. Model delimiters and
        // the score remarks of a previous docking run are not part of the
        // molecule; the writer emits fresh ones.
        if !line.starts_with("MODEL") && !line.starts_with("ENDMDL") && !is_score_remark(&line) {
            molecule.pdbqt_records.push(line.clone());
        }

        // Parse based on record type
        if line.starts_with("ATOM") || line.starts_with("HETATM") {
            // Parse atom record
            let atom = parse_pdbqt_atom(&line, line_number)?;
            let atom_idx = molecule.add_atom(atom);
            fragment_of.push(*fragment_stack.last().unwrap_or(&0));

            // Add to current branch if any
            if !branch_stack.is_empty() {
                current_branch_atoms.push(atom_idx);
            }
        } else if line.starts_with("MODEL") {
            // Start a new model
            if model_number > 1 {
                // Only use the first model for now
                break;
            }
            model_number += 1;
        } else if line.starts_with("ENDMDL") {
            // End of a model
            break;
        } else if line.starts_with("BRANCH") {
            // Parse BRANCH record
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() != 3 {
                return Err(IoError::Parse {
                    line: line_number,
                    message: format!("Invalid BRANCH record: {}", line),
                });
            }

            let from_atom: usize = parts[1].parse().map_err(|_| IoError::Parse {
                line: line_number,
                message: format!("Invalid BRANCH from atom: {}", parts[1]),
            })?;

            let to_atom: usize = parts[2].parse().map_err(|_| IoError::Parse {
                line: line_number,
                message: format!("Invalid BRANCH to atom: {}", parts[2]),
            })?;

            // Convert 1-based indices from PDBQT to 0-based indices
            let from_atom = from_atom - 1;
            let to_atom = to_atom - 1;

            // Push current branch atoms onto stack (if any)
            if !current_branch_atoms.is_empty() {
                branch_atoms.push(current_branch_atoms.clone());
            }

            // Start a new branch - store atom indices for later bond creation
            branch_stack.push((from_atom, to_atom));
            current_branch_atoms = Vec::new();
            fragment_stack.push(fragment_count);
            fragment_count += 1;
        } else if line.starts_with("ENDBRANCH") {
            // Parse ENDBRANCH record
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() != 3 {
                return Err(IoError::Parse {
                    line: line_number,
                    message: format!("Invalid ENDBRANCH record: {}", line),
                });
            }
            if fragment_stack.len() > 1 {
                fragment_stack.pop();
            }

            if let Some((from_atom, to_atom)) = branch_stack.pop() {
                // Now both atoms should exist - create the rotatable bond
                if from_atom < molecule.atoms.len() && to_atom < molecule.atoms.len() {
                    // Add the rotatable bond
                    let bond_idx = molecule.bonds.len();
                    molecule.bonds.push(Bond {
                        atom1_idx: from_atom,
                        atom2_idx: to_atom,
                        rotatable: true,
                    });

                    // Create a torsion for this rotatable bond
                    let torsion = Torsion {
                        bond_idx,
                        moving_atoms: current_branch_atoms.clone(),
                        angle: 0.0, // Default angle
                    };
                    molecule.torsions.push(torsion);
                }

                // Restore previous branch's atoms and add current branch atoms to it
                let branch_atom_copy = current_branch_atoms.clone();
                if !branch_atoms.is_empty() {
                    current_branch_atoms = branch_atoms.pop().unwrap();
                    // Add atoms from this branch to parent branch
                    current_branch_atoms.extend(branch_atom_copy);
                } else {
                    current_branch_atoms = Vec::new();
                }
            }
        } else if line.starts_with("TORSDOF") {
            // Number of torsional degrees of freedom
            // This is just informational, but we could use it to validate
        } else if line.starts_with("REMARK") {
            // Ignore remarks for now
        } else if line.starts_with("ROOT") || line.starts_with("ENDROOT") {
            // These mark the start and end of the rigid part of the ligand
        } else {
            // Unrecognized record type - ignore
        }
    }

    // Always infer connectivity. BRANCH records only describe the rotatable
    // bonds (6 for a typical drug-like ligand), which is nowhere near the full
    // bond graph that atom typing needs.
    infer_bonds(&mut molecule, &fragment_of);

    // Hydrophobicity and donor status depend on connectivity, so they can only
    // be assigned once bonds exist.
    molecule.assign_xs_types();

    Ok(molecule)
}

/// Parse an atom record from a PDBQT file
fn parse_pdbqt_atom(line: &str, line_number: usize) -> Result<Atom, IoError> {
    if line.len() < 78 {
        return Err(IoError::Parse {
            line: line_number,
            message: format!("Line too short for atom record: {}", line),
        });
    }

    // Parse atom serial number
    let serial = line[6..11]
        .trim()
        .parse::<u32>()
        .map_err(|_| IoError::Parse {
            line: line_number,
            message: format!("Invalid atom serial number: {}", &line[6..11]),
        })?;

    // Parse atom name
    let name = line[12..16].trim().to_string();

    // Parse residue name
    let residue_name = line[17..20].trim().to_string();

    // Parse chain ID
    let chain_id = line[21..22].trim().chars().next().unwrap_or('A');

    // Parse residue number
    let residue_num = line[22..26]
        .trim()
        .parse::<i32>()
        .map_err(|_| IoError::Parse {
            line: line_number,
            message: format!("Invalid residue number: {}", &line[22..26]),
        })?;

    // Parse coordinates
    let x = line[30..38]
        .trim()
        .parse::<f64>()
        .map_err(|_| IoError::Parse {
            line: line_number,
            message: format!("Invalid x coordinate: {}", &line[30..38]),
        })?;

    let y = line[38..46]
        .trim()
        .parse::<f64>()
        .map_err(|_| IoError::Parse {
            line: line_number,
            message: format!("Invalid y coordinate: {}", &line[38..46]),
        })?;

    let z = line[46..54]
        .trim()
        .parse::<f64>()
        .map_err(|_| IoError::Parse {
            line: line_number,
            message: format!("Invalid z coordinate: {}", &line[46..54]),
        })?;

    // Parse partial charge (optional)
    let charge = if line.len() >= 76 {
        line[68..76].trim().parse::<f64>().unwrap_or(0.0)
    } else {
        0.0
    };

    // Parse atom type: column 78 to the end of the line, as Vina reads it,
    // since some types (`CG0`) run past the two-column field.
    let atom_type_str = line.get(77..).unwrap_or("").trim();

    let atom_type = AtomType::from_pdbqt_string(atom_type_str);

    // Create atom
    let atom = Atom::new(
        atom_type,
        Vector3::new(x, y, z),
        name,
        serial,
        residue_name,
        residue_num,
        chain_id,
        charge,
    );

    Ok(atom)
}

/// Infer bonds based on distances between atoms.
///
/// As in Vina, only atoms of the same rigid fragment (`fragment_of`) can
/// bond by distance; bonds between fragments are the BRANCH bonds. This
/// keeps an opened macrocycle open: Meeko puts each closure dummy right on
/// top of a carbon in another fragment.
fn infer_bonds(molecule: &mut Molecule, fragment_of: &[usize]) {
    // Heavy-atom covalent bonds run to about 1.9 A (S-S is the long case at
    // ~2.05); bonds to hydrogen are near 1.1. Using one flat threshold for both
    // would pull in 1-3 neighbours, which sit around 2.4 A.
    const HEAVY_BOND_MAX: f64 = 2.1;
    const HYDROGEN_BOND_MAX: f64 = 1.3;

    // Rotatable bonds come from the BRANCH records and are already present;
    // record them so inference does not duplicate them.
    let mut existing: HashSet<(usize, usize)> = molecule
        .bonds
        .iter()
        .map(|b| ordered_pair(b.atom1_idx, b.atom2_idx))
        .collect();

    for i in 0..molecule.atoms.len() {
        for j in (i + 1)..molecule.atoms.len() {
            if fragment_of.get(i) != fragment_of.get(j) {
                continue;
            }
            let atom1 = &molecule.atoms[i];
            let atom2 = &molecule.atoms[j];

            let involves_hydrogen = atom1.atom_type.is_hydrogen() || atom2.atom_type.is_hydrogen();
            let threshold = if involves_hydrogen {
                HYDROGEN_BOND_MAX
            } else {
                HEAVY_BOND_MAX
            };

            if atom1.distance(atom2) >= threshold {
                continue;
            }
            if !existing.insert(ordered_pair(i, j)) {
                continue;
            }

            // Inferred bonds describe connectivity only. Rotatability is a
            // property of the torsion tree, which the BRANCH records define;
            // guessing it from distance alone would mark nearly every bond
            // rotatable.
            let _ = molecule.add_bond(i, j, false);
        }
    }
}

fn ordered_pair(i: usize, j: usize) -> (usize, usize) {
    if i <= j {
        (i, j)
    } else {
        (j, i)
    }
}

/// Whether `line` is one of the per-pose score remarks Vina writes.
fn is_score_remark(line: &str) -> bool {
    [
        "REMARK VINA RESULT",
        "REMARK INTER",
        "REMARK INTRA",
        "REMARK UNBOUND",
    ]
    .iter()
    .any(|prefix| line.starts_with(prefix))
}

/// Write a molecule to a PDBQT file.
///
/// A molecule read from PDBQT is written back record for record, with only
/// its coordinates changed (see [`Molecule::pdbqt_records`]).
pub fn write_pdbqt<P: AsRef<Path>>(molecule: &Molecule, path: P) -> Result<(), IoError> {
    let mut file = File::create(path)?;
    write_molecule_records(&mut file, molecule)
}

/// Write docking results as a multi-model PDBQT file in Vina's layout.
///
/// Each model carries Vina's score remarks (`VINA RESULT` with the affinity
/// and RMSD from the best mode, then the inter/intra breakdown when the
/// result has it) followed by the ligand's own records, so Vina, Meeko and
/// rustdock all read it back.
pub fn write_docking_results<P: AsRef<Path>>(
    results: &[crate::optimization::DockingResult],
    path: P,
) -> Result<(), IoError> {
    let mut file = File::create(path)?;

    for (i, result) in results.iter().enumerate() {
        writeln!(file, "MODEL {}", i + 1)?;
        let rmsd = result.rmsd.unwrap_or(0.0);
        writeln!(
            file,
            "REMARK VINA RESULT: {:9.3} {:10.3} {:10.3}",
            result.energy, rmsd, rmsd
        )?;
        let component = |name: &str| {
            result
                .energy_components
                .iter()
                .find(|(n, _)| n == name)
                .map(|(_, v)| *v)
        };
        if let (Some(inter), Some(intra)) =
            (component("Intermolecular"), component("Intramolecular"))
        {
            writeln!(
                file,
                "{:<21}{:>16.3}",
                "REMARK INTER + INTRA:",
                inter + intra
            )?;
            writeln!(file, "{:<21}{:>16.3}", "REMARK INTER:", inter)?;
            writeln!(file, "{:<21}{:>16.3}", "REMARK INTRA:", intra)?;
            if let Some(unbound) = component("Unbound") {
                writeln!(file, "{:<21}{:>16.3}", "REMARK UNBOUND:", unbound)?;
            }
        }
        write_molecule_records(&mut file, &result.molecule)?;
        writeln!(file, "ENDMDL")?;
    }

    Ok(())
}

/// Write a molecule's records: its original PDBQT records with coordinates
/// substituted when it has them, otherwise plain atom records.
fn write_molecule_records<W: Write>(out: &mut W, molecule: &Molecule) -> Result<(), IoError> {
    let is_atom = |l: &str| l.starts_with("ATOM") || l.starts_with("HETATM");
    let records = &molecule.pdbqt_records;
    let atom_records = records.iter().filter(|l| is_atom(l)).count();

    if !records.is_empty() && atom_records == molecule.atoms.len() {
        let mut atoms = molecule.atoms.iter();
        for line in records {
            if is_atom(line) {
                let atom = atoms.next().expect("counted above");
                writeln!(out, "{}", with_coordinates(line, &atom.coordinates))?;
            } else {
                writeln!(out, "{line}")?;
            }
        }
        return Ok(());
    }

    for (i, atom) in molecule.atoms.iter().enumerate() {
        let record = if atom.name.starts_with("HETATM") {
            "HETATM"
        } else {
            "ATOM  "
        };
        // PDB convention: names shorter than four characters start one
        // column in, so the element lines up.
        let name = if atom.name.len() < 4 {
            format!(" {:<3}", atom.name)
        } else {
            format!("{:<4}", atom.name)
        };
        writeln!(
            out,
            "{}{:5} {} {:>3} {:1}{:4}    {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}    {:6.3} {:<2}",
            record,
            (i + 1) % 100_000,
            name,
            atom.residue_name,
            atom.chain_id,
            atom.residue_num,
            atom.coordinates.x,
            atom.coordinates.y,
            atom.coordinates.z,
            1.0,
            0.0,
            atom.charge,
            atom.atom_type.to_pdbqt_string()
        )?;
    }
    Ok(())
}

/// An atom record with its coordinate columns (31-54) replaced.
fn with_coordinates(line: &str, p: &Vector3<f64>) -> String {
    let head = line.get(..30).unwrap_or(line);
    let tail = line.get(54..).unwrap_or("");
    format!("{head:<30}{:8.3}{:8.3}{:8.3}{tail}", p.x, p.y, p.z)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("test_data")
            .join(name)
    }

    /// A Meeko ligand written back unchanged must be byte-identical: its
    /// SMILES remarks, torsion tree and atom columns all survive.
    #[test]
    fn test_meeko_ligand_round_trips_exactly() {
        let path = fixture("ligand_meeko.pdbqt");
        let ligand = parse_pdbqt(&path).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.pdbqt");
        write_pdbqt(&ligand, &out).unwrap();
        assert_eq!(
            std::fs::read_to_string(&out).unwrap().trim_end(),
            std::fs::read_to_string(&path).unwrap().trim_end()
        );
    }

    /// Docking output keeps every model's ligand records, gains Vina's score
    /// remarks, and reads back to the first pose; re-reading does not carry
    /// the old score remarks into the molecule.
    #[test]
    fn test_docking_results_read_back() {
        let ligand = parse_pdbqt(fixture("ligand_meeko.pdbqt")).unwrap();
        let mut moved = ligand.clone();
        moved.translate(&Vector3::new(1.0, -2.0, 0.5));
        let result = |molecule: Molecule, energy: f64| crate::optimization::DockingResult {
            energy,
            molecule,
            energy_components: vec![
                ("Intermolecular".to_string(), -9.0),
                ("Intramolecular".to_string(), -0.5),
                ("Unbound".to_string(), -0.5),
            ],
            rmsd: None,
        };
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("docked.pdbqt");
        write_docking_results(
            &[result(moved.clone(), -7.25), result(ligand.clone(), -6.0)],
            &out,
        )
        .unwrap();

        let text = std::fs::read_to_string(&out).unwrap();
        assert!(text.starts_with("MODEL 1\nREMARK VINA RESULT:    -7.250      0.000      0.000\n"));
        assert!(text.contains("REMARK INTER + INTRA:          -9.500\n"));
        assert_eq!(text.matches("REMARK SMILES Cc1").count(), 2);
        assert_eq!(text.matches("TORSDOF").count(), 2);

        let back = parse_pdbqt(&out).unwrap();
        assert_eq!(back.atoms.len(), moved.atoms.len());
        assert_eq!(back.torsions.len(), moved.torsions.len());
        for (a, b) in back.atoms.iter().zip(&moved.atoms) {
            assert!((a.coordinates - b.coordinates).norm() < 1e-3);
        }
        assert!(back.pdbqt_records.iter().all(|l| !is_score_remark(l)));
        assert_eq!(back.pdbqt_records.len(), ligand.pdbqt_records.len());
    }

    /// A six-membered ring opened as Meeko writes it: closure carbons 1 and
    /// 7 in different rigid fragments, each carrying a dummy (2, 8) that
    /// sits exactly on the other.
    const OPENED_RING: &str = "\
ROOT
ATOM      1  C   UNL     1       1.500   0.000   0.000  1.00  0.00     0.000 CG0
ATOM      2  G   UNL     1       0.750  -1.300   0.000  1.00  0.00     0.000 G0
ATOM      3  C   UNL     1       0.750   1.300   0.000  1.00  0.00     0.000 C
ENDROOT
BRANCH   3   4
ATOM      4  C   UNL     1      -0.750   1.300   0.000  1.00  0.00     0.000 C
BRANCH   4   5
ATOM      5  C   UNL     1      -1.500   0.000   0.000  1.00  0.00     0.000 C
BRANCH   5   6
ATOM      6  C   UNL     1      -0.750  -1.300   0.000  1.00  0.00     0.000 C
ATOM      7  C   UNL     1       0.750  -1.300   0.000  1.00  0.00     0.000 CG0
ATOM      8  G   UNL     1       1.500   0.000   0.000  1.00  0.00     0.000 G0
ENDBRANCH   5   6
ENDBRANCH   4   5
ENDBRANCH   3   4
TORSDOF 3
";

    /// Closure types parse in full, and bonds form only within a rigid
    /// fragment, so the ring stays open and each dummy is glued, not
    /// bonded, to its partner.
    #[test]
    fn test_opened_macrocycle() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("ring.pdbqt");
        std::fs::write(&path, OPENED_RING).unwrap();
        let ring = parse_pdbqt(&path).unwrap();

        assert_eq!(ring.atoms[0].atom_type, AtomType::CarbonClosure(0));
        assert_eq!(ring.atoms[1].atom_type, AtomType::ClosureDummy(0));
        let mut bonds: Vec<_> = ring
            .bonds
            .iter()
            .map(|b| ordered_pair(b.atom1_idx, b.atom2_idx))
            .collect();
        bonds.sort();
        assert_eq!(
            bonds,
            [(0, 1), (0, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)]
        );
        // The carbons next to a closure carbon are polar, as Vina has it.
        assert!(!ring.atoms[2].is_hydrophobic && !ring.atoms[5].is_hydrophobic);
        assert!(ring.atoms[3].is_hydrophobic);

        assert_eq!(ring.glue_pairs(), [(0, 7), (1, 6)]);
        // Closed, every heavy pair is within three bonds of the closure, and
        // dummies pair with nothing else.
        assert!(ring.intramolecular_pairs().is_empty());
    }

    #[test]
    fn test_negative_residue_number() {
        // PDB allows negative residue numbers (e.g. tags before residue 1),
        // and Open Babel writes them through to receptor PDBQTs.
        let line =
            "ATOM      1  N   PHE A  -1      21.700  33.210  12.363  0.00  0.00    +0.000 NA";
        let atom = parse_pdbqt_atom(line, 1).unwrap();
        assert_eq!(atom.residue_num, -1);
        assert_eq!(atom.residue_name, "PHE");
    }
}
