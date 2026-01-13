//! Input/output functionality for molecular docking

use nalgebra::Vector3;
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

        // Parse based on record type
        if line.starts_with("ATOM") || line.starts_with("HETATM") {
            // Parse atom record
            let atom = parse_pdbqt_atom(&line, line_number)?;
            let atom_idx = molecule.add_atom(atom);

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
        } else if line.starts_with("ENDBRANCH") {
            // Parse ENDBRANCH record
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() != 3 {
                return Err(IoError::Parse {
                    line: line_number,
                    message: format!("Invalid ENDBRANCH record: {}", line),
                });
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

    // If no bonds were created, try to infer them based on distances
    if molecule.bonds.is_empty() {
        infer_bonds(&mut molecule);
    }

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
        .parse::<u32>()
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

    // Parse atom type
    let atom_type_str = if line.len() >= 79 {
        line[77..79].trim()
    } else {
        ""
    };

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

/// Infer bonds based on distances between atoms
fn infer_bonds(molecule: &mut Molecule) {
    // Typical bond length threshold in Angstroms
    const BOND_THRESHOLD: f64 = 2.0;

    for i in 0..molecule.atoms.len() {
        for j in (i + 1)..molecule.atoms.len() {
            let atom1 = &molecule.atoms[i];
            let atom2 = &molecule.atoms[j];

            let distance = atom1.distance(atom2);

            // If atoms are close enough, add a bond
            if distance < BOND_THRESHOLD {
                // Determine if the bond is rotatable
                // This is a simplification - in reality, this would be more complex
                let rotatable =
                    !matches!(atom1.atom_type, AtomType::Hydrogen | AtomType::HydrogenD)
                        && !matches!(atom2.atom_type, AtomType::Hydrogen | AtomType::HydrogenD);

                let _ = molecule.add_bond(i, j, rotatable);
            }
        }
    }
}

/// Write a molecule to a PDBQT file
pub fn write_pdbqt<P: AsRef<Path>>(molecule: &Molecule, path: P) -> Result<(), IoError> {
    let mut file = File::create(path)?;

    // Write header
    writeln!(file, "REMARK PDBQT file generated by rustdock-vina")?;

    // Write atom records
    for (i, atom) in molecule.atoms.iter().enumerate() {
        let record_type = if atom.name.starts_with("HETATM") {
            "HETATM"
        } else {
            "ATOM  "
        };

        writeln!(
            file,
            "{}{:5} {:4} {:3} {:1}{:4}    {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}     {:8.3} {}",
            record_type,
            i + 1, // 1-based index
            atom.name,
            atom.residue_name,
            atom.chain_id,
            atom.residue_num,
            atom.coordinates.x,
            atom.coordinates.y,
            atom.coordinates.z,
            1.0, // Occupancy
            0.0, // Temperature factor
            atom.charge,
            atom.atom_type.to_pdbqt_string()
        )?;
    }

    // Write trailer
    writeln!(file, "END")?;

    Ok(())
}

/// Write docking results to a PDBQT file
pub fn write_docking_results<P: AsRef<Path>>(
    results: &[crate::optimization::DockingResult],
    path: P,
) -> Result<(), IoError> {
    let mut file = File::create(path)?;

    // Write header
    writeln!(file, "REMARK PDBQT file generated by rustdock-vina")?;
    writeln!(file, "REMARK Contains {} docked models", results.len())?;

    for (i, result) in results.iter().enumerate() {
        // Write model header
        writeln!(file, "MODEL {:>4}", i + 1)?;
        writeln!(
            file,
            "REMARK VINA RESULT:    {:6.3} kcal/mol",
            result.energy
        )?;

        if let Some(rmsd) = result.rmsd {
            writeln!(file, "REMARK RMSD FROM BEST MODE: {:6.3} Å", rmsd)?;
        }

        // Write atom records
        for (j, atom) in result.molecule.atoms.iter().enumerate() {
            let record_type = if atom.name.starts_with("HETATM") {
                "HETATM"
            } else {
                "ATOM  "
            };

            writeln!(
                file,
                "{}{:5} {:4} {:3} {:1}{:4}    {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}     {:8.3} {}",
                record_type,
                j + 1, // 1-based index
                atom.name,
                atom.residue_name,
                atom.chain_id,
                atom.residue_num,
                atom.coordinates.x,
                atom.coordinates.y,
                atom.coordinates.z,
                1.0, // Occupancy
                0.0, // Temperature factor
                atom.charge,
                atom.atom_type.to_pdbqt_string()
            )?;
        }

        // Write model trailer
        writeln!(file, "ENDMDL")?;
    }

    // Write file trailer
    writeln!(file, "END")?;

    Ok(())
}
