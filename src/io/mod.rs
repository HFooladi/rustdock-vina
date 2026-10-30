//! Input/output functionality for molecular docking

use nalgebra::Vector3;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use thiserror::Error;

use crate::atom::{Atom, AtomType};
use crate::molecule::{Bond, FlexResidues, Molecule, Torsion};

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

    let mut tree = TreeBuilder::default();
    let mut fragments = Fragments::default();

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
            let serial = atom.serial;
            let atom_idx = molecule.add_atom(atom);
            tree.atom(atom_idx, serial, &mut fragments);
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
            tree.branch(&line, line_number, &mut fragments)?;
        } else if line.starts_with("ENDBRANCH") {
            tree.end_branch(&line, line_number, &mut molecule, &mut fragments)?;
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
    infer_bonds(&mut molecule, |a, b| fragments.rigid(a, b));

    // Hydrophobicity and donor status depend on connectivity, so they can only
    // be assigned once bonds exist.
    molecule.assign_xs_types();

    Ok(molecule)
}

/// Parse a flexible residue PDBQT (Vina's `--flex`) for the rigid receptor
/// `rigid`.
///
/// Each `BEGIN_RES`…`END_RES` block has a ROOT and BRANCHes, each one a
/// torsion. The ROOT atoms stay fixed, and so does the first atom of each
/// branch off the ROOT, which lies on its rotation axis. The fixed atoms
/// are bonded to the rigid receptor's atoms as well as the residue's own,
/// as Vina does, so a CA whose N is in the rigid file is typed as bonded to
/// a heteroatom.
pub fn parse_flex_pdbqt<P: AsRef<Path>>(
    path: P,
    rigid: &Molecule,
) -> Result<FlexResidues, IoError> {
    let file = File::open(path.as_ref())?;
    let mut molecule = Molecule::new(
        path.as_ref()
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("flex"),
    );
    let mut residue_of = Vec::new();
    let mut fragments = Fragments::default();
    let mut tree = TreeBuilder::default();
    let mut residue: Option<usize> = None;
    let mut residues = 0;

    for (line_number, line) in BufReader::new(file).lines().enumerate() {
        let line = line?;
        let line_number = line_number + 1;
        if line.trim().is_empty() {
            continue;
        }
        if !is_score_remark(&line) {
            molecule.pdbqt_records.push(line.clone());
        }
        let outside = |what: &str| IoError::Parse {
            line: line_number,
            message: format!("{what} outside a BEGIN_RES/END_RES block: {line}"),
        };

        if line.starts_with("BEGIN_RES") {
            if residue.is_some() {
                return Err(IoError::Parse {
                    line: line_number,
                    message: "BEGIN_RES before the previous END_RES".to_string(),
                });
            }
            residue = Some(residues);
            residues += 1;
            // Serials and open branches are per residue; fragment numbers
            // are not, and every residue's ROOT is in fragment 0, as all of
            // Vina's fixed flex atoms are rigid relative to each other.
            tree = TreeBuilder::default();
        } else if line.starts_with("END_RES") {
            if residue.take().is_none() || !tree.open.is_empty() {
                return Err(IoError::Parse {
                    line: line_number,
                    message: format!("unbalanced END_RES: {line}"),
                });
            }
        } else if line.starts_with("ATOM") || line.starts_with("HETATM") {
            let index = residue.ok_or_else(|| outside("atom"))?;
            let atom = parse_pdbqt_atom(&line, line_number)?;
            let serial = atom.serial;
            let atom_idx = molecule.add_atom(atom);
            residue_of.push(index);
            tree.atom(atom_idx, serial, &mut fragments);
        } else if line.starts_with("ENDBRANCH") {
            residue.ok_or_else(|| outside("ENDBRANCH"))?;
            tree.end_branch(&line, line_number, &mut molecule, &mut fragments)?;
        } else if line.starts_with("BRANCH") {
            residue.ok_or_else(|| outside("BRANCH"))?;
            tree.branch(&line, line_number, &mut fragments)?;
        } else if line.starts_with("MODEL") || line.starts_with("ENDMDL") {
            return Err(IoError::InvalidFormat(format!(
                "{}: a flex residue file has one model; split it first",
                path.as_ref().display()
            )));
        }
        // ROOT/ENDROOT, REMARK and anything else carry nothing to model.
    }
    if residue.is_some() {
        return Err(IoError::InvalidFormat(format!(
            "{}: missing END_RES",
            path.as_ref().display()
        )));
    }
    molecule.is_rigid = molecule.torsions.is_empty();

    // Vina's rigid bodies: a branch's first atom lies on its rotation axis,
    // so it moves with the parent's body, not its own. The ROOT's body,
    // with the first atom of every branch off it, never moves.
    let body_of: Vec<usize> = (0..molecule.atoms.len())
        .map(|i| {
            let f = fragments.of[i];
            match fragments.parent[f] {
                Some(parent) if fragments.first[f] == Some(i) => fragments.of[parent],
                _ => f,
            }
        })
        .collect();
    let mut flex = FlexResidues {
        fixed: body_of.iter().map(|&b| b == 0).collect(),
        molecule,
        residue_of,
        body_of,
        axis: fragments.parent.into_iter().zip(fragments.first).collect(),
    };
    type_flex_atoms(&mut flex, rigid);
    Ok(flex)
}

/// Infer the flex residues' bonds and assign their XS types.
///
/// Vina infers bonds over the rigid and flex atoms together, between atoms
/// that cannot move relative to each other. The fixed flex atoms then also
/// bond to the rigid receptor's atoms, which can make a fixed carbon polar.
/// The rigid atoms near a fixed atom join the bond inference; only the
/// bonds among flex atoms are kept.
fn type_flex_atoms(flex: &mut FlexResidues, rigid: &Molecule) {
    // Far enough for a bond to a fixed atom, and for an atom between them.
    const REACH: f64 = 4.0;
    let n = flex.molecule.atoms.len();
    let fixed: Vec<_> = flex
        .molecule
        .atoms
        .iter()
        .zip(&flex.fixed)
        .filter(|(_, &f)| f)
        .map(|(a, _)| a.coordinates)
        .collect();
    let mut combined = flex.molecule.clone();
    for atom in &rigid.atoms {
        if fixed
            .iter()
            .any(|r| (r - atom.coordinates).norm_squared() < REACH * REACH)
        {
            combined.add_atom(atom.clone());
        }
    }
    // The rigid atoms cannot move relative to each other or to the fixed
    // flex atoms.
    let still = |i: usize| i >= n || flex.fixed[i];
    infer_bonds(&mut combined, |a, b| {
        if a < n && b < n {
            !flex.variable(a, b)
        } else {
            still(a) && still(b)
        }
    });
    combined.assign_xs_types();

    combined.atoms.truncate(n);
    combined
        .bonds
        .retain(|b| b.atom1_idx < n && b.atom2_idx < n);
    flex.molecule.atoms = combined.atoms;
    flex.molecule.bonds = combined.bonds;
}

/// The rigid fragments of a torsion tree, which decide which atoms can bond
/// and which pairs can move relative to each other.
struct Fragments {
    /// Fragment of each atom: the ROOT is 0, each BRANCH opens a new one.
    of: Vec<usize>,
    /// The atom each fragment hangs from (the first atom of its BRANCH
    /// record).
    parent: Vec<Option<usize>>,
    /// Each fragment's first atom, bonded to its parent (the second atom of
    /// its BRANCH record), once the branch is closed.
    first: Vec<Option<usize>>,
}

impl Fragments {
    /// Whether atoms `a` and `b` cannot move relative to each other: they
    /// are in the same fragment, or one is the atom the other's fragment
    /// hangs from.
    fn rigid(&self, a: usize, b: usize) -> bool {
        match (self.of.get(a), self.of.get(b)) {
            (Some(&fa), Some(&fb)) => {
                fa == fb
                    || self.parent.get(fb) == Some(&Some(a))
                    || self.parent.get(fa) == Some(&Some(b))
            }
            _ => false,
        }
    }
}

impl Default for Fragments {
    fn default() -> Self {
        Self {
            of: Vec::new(),
            parent: vec![None],
            first: vec![None],
        }
    }
}

/// Builds a ligand's or flexible residue's torsions from its BRANCH and
/// ENDBRANCH records.
///
/// BRANCH records name atoms by serial number, which Vina looks up among the
/// atoms read so far, so a flex residue file keeps the receptor's serials.
struct TreeBuilder {
    /// Open branches: the atom each hangs from, and the serial of its first
    /// atom (not yet read when the BRANCH record is).
    open: Vec<(usize, u32)>,
    /// The moving atoms of each open branch, innermost last.
    moving: Vec<Vec<usize>>,
    /// The fragment of each open branch, with the ROOT's at the bottom.
    fragment_stack: Vec<usize>,
    /// Atom index of each serial in the current ligand or residue.
    serials: HashMap<u32, usize>,
}

impl Default for TreeBuilder {
    fn default() -> Self {
        Self {
            open: Vec::new(),
            moving: Vec::new(),
            fragment_stack: vec![0],
            serials: HashMap::new(),
        }
    }
}

impl TreeBuilder {
    /// Record the atom just added at `index`.
    fn atom(&mut self, index: usize, serial: u32, fragments: &mut Fragments) {
        fragments.of.push(*self.fragment_stack.last().unwrap_or(&0));
        self.serials.insert(serial, index);
        // An atom moves with every enclosing branch; ENDBRANCH passes a
        // branch's atoms on to its parent.
        if let Some(moving) = self.moving.last_mut() {
            moving.push(index);
        }
    }

    fn branch(
        &mut self,
        line: &str,
        line_number: usize,
        fragments: &mut Fragments,
    ) -> Result<(), IoError> {
        let (from, to) = branch_serials(line, line_number)?;
        let from = self.lookup(from, line, line_number)?;
        self.open.push((from, to));
        self.moving.push(Vec::new());
        self.fragment_stack.push(fragments.parent.len());
        fragments.parent.push(Some(from));
        fragments.first.push(None);
        Ok(())
    }

    /// Close the innermost branch: add its rotatable bond and torsion.
    fn end_branch(
        &mut self,
        line: &str,
        line_number: usize,
        molecule: &mut Molecule,
        fragments: &mut Fragments,
    ) -> Result<(), IoError> {
        branch_serials(line, line_number)?;
        let (Some((from, to)), Some(moving)) = (self.open.pop(), self.moving.pop()) else {
            return Err(IoError::Parse {
                line: line_number,
                message: format!("ENDBRANCH without a BRANCH: {line}"),
            });
        };
        let fragment = self.fragment_stack.pop().unwrap_or(0);
        let to = self.lookup(to, line, line_number)?;
        fragments.first[fragment] = Some(to);

        let bond_idx = molecule.bonds.len();
        molecule.bonds.push(Bond {
            atom1_idx: from,
            atom2_idx: to,
            rotatable: true,
        });
        if let Some(parent) = self.moving.last_mut() {
            parent.extend_from_slice(&moving);
        }
        molecule.torsions.push(Torsion {
            bond_idx,
            moving_atoms: moving,
            angle: 0.0,
        });
        Ok(())
    }

    fn lookup(&self, serial: u32, line: &str, line_number: usize) -> Result<usize, IoError> {
        self.serials
            .get(&serial)
            .copied()
            .ok_or_else(|| IoError::Parse {
                line: line_number,
                message: format!("atom {serial} of this branch has not been read: {line}"),
            })
    }
}

/// The two serial numbers of a BRANCH or ENDBRANCH record.
fn branch_serials(line: &str, line_number: usize) -> Result<(u32, u32), IoError> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    let number = |s: &str| s.parse::<u32>().ok();
    match parts.as_slice() {
        [_, a, b] => number(a).zip(number(b)),
        _ => None,
    }
    .ok_or_else(|| IoError::Parse {
        line: line_number,
        message: format!("Invalid branch record: {line}"),
    })
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

/// Infer bonds from distances, as Vina's `model::assign_bonds` does.
///
/// Two atoms can only bond if they cannot move relative to each other
/// (`rigid`); in a ligand, they are in the same rigid fragment, or one is
/// the atom the other's fragment hangs from. This keeps an opened
/// macrocycle open, since Meeko puts each closure dummy right on top of a
/// carbon in another fragment. They bond if closer than 1.1 times the sum
/// of their covalent radii, unless a third such atom is closer to both of
/// them than they are to each other. The bonds between fragments are the
/// BRANCH bonds, which are already present.
fn infer_bonds(molecule: &mut Molecule, rigid: impl Fn(usize, usize) -> bool) {
    const ALLOWANCE: f64 = 1.1;
    // The largest covalent radius; Vina looks this far for atoms in between.
    const MAX_COVALENT: f64 = 1.75;

    let mut existing: HashSet<(usize, usize)> = molecule
        .bonds
        .iter()
        .map(|b| ordered_pair(b.atom1_idx, b.atom2_idx))
        .collect();

    // Bond candidates, and each atom's rigidly attached neighbours that could
    // stand between it and a bonding partner.
    let n = molecule.atoms.len();
    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut candidates = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if !rigid(i, j) {
                continue;
            }
            let (a, b) = (&molecule.atoms[i], &molecule.atoms[j]);
            let r2 = (a.coordinates - b.coordinates).norm_squared();
            let (ra, rb) = (a.atom_type.covalent_radius(), b.atom_type.covalent_radius());
            let reach = |r: f64| (ALLOWANCE * (r + MAX_COVALENT)).powi(2);
            if r2 < reach(ra) {
                neighbours[i].push(j);
            }
            if r2 < reach(rb) {
                neighbours[j].push(i);
            }
            if r2 < (ALLOWANCE * (ra + rb)).powi(2) {
                candidates.push((i, j, r2));
            }
        }
    }

    for (i, j, r2) in candidates {
        let between = neighbours[i].iter().any(|&c| {
            c != j
                && rigid(j, c)
                && (molecule.atoms[i].coordinates - molecule.atoms[c].coordinates).norm_squared()
                    < r2
                && (molecule.atoms[j].coordinates - molecule.atoms[c].coordinates).norm_squared()
                    < r2
        });
        if between || !existing.insert(ordered_pair(i, j)) {
            continue;
        }
        // Inferred bonds describe connectivity only. Rotatability is a
        // property of the torsion tree, which the BRANCH records define;
        // guessing it from distance alone would mark nearly every bond
        // rotatable.
        let _ = molecule.add_bond(i, j, false);
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

    /// The N-methyl of 2wnc's ligand, whose amine hydrogen sits closer to both
    /// the N and the methyl carbon than they are to each other. Vina's
    /// `atom_exists_between` then leaves the N-C bond out, so the carbon is
    /// hydrophobic.
    #[test]
    fn test_atom_between_blocks_bond() {
        let text = "\
ROOT
ATOM      1  C   UNL     1     -70.593 -47.015 -74.539  1.00  0.00     0.187 C
ATOM      2  N   UNL     1     -71.126 -45.907 -73.741  1.00  0.00    -0.332 N
ATOM      3  C   UNL     1     -70.670 -44.553 -73.842  1.00  0.00     0.333 C
ATOM      4  H   UNL     1     -71.725 -45.282 -74.420  1.00  0.00     0.348 HD
ENDROOT
TORSDOF 0
";
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("amine.pdbqt");
        std::fs::write(&path, text).unwrap();
        let amine = parse_pdbqt(&path).unwrap();

        let mut bonds: Vec<_> = amine
            .bonds
            .iter()
            .map(|b| ordered_pair(b.atom1_idx, b.atom2_idx))
            .collect();
        bonds.sort();
        assert_eq!(bonds, [(0, 1), (1, 3)]);
        assert!(!amine.atoms[0].is_hydrophobic);
        assert!(amine.atoms[2].is_hydrophobic);
        assert!(amine.atoms[1].is_hbond_donor);
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

    /// Two residues as Meeko writes them, with the serials renumbered so
    /// that BRANCH records cannot be read as file positions.
    const FLEX: &str = "\
BEGIN_RES TYR A 79
REMARK INDEX MAP 3 1 15 2 18 3 19 4 20 5 22 6 23 7 26 8 28 9 29 10
ROOT
ATOM    101  CA  TYR A  79      25.291   6.625  63.328  1.00  0.00     0.175 C 
ENDROOT
BRANCH 101 102
ATOM    102  CB  TYR A  79      25.190   5.178  62.860  1.00  0.00     0.070 C 
BRANCH 102 103
ATOM    103  CG  TYR A  79      25.993   4.196  63.686  1.00  0.00    -0.045 A 
ATOM    104 CD1  TYR A  79      27.363   4.036  63.481  1.00  0.00     0.007 A 
ATOM    105 CD2  TYR A  79      25.381   3.380  64.599  1.00  0.00     0.007 A 
ATOM    106 CE1  TYR A  79      28.091   3.119  64.198  1.00  0.00     0.046 A 
ATOM    107 CE2  TYR A  79      26.115   2.443  65.353  1.00  0.00     0.046 A 
ATOM    108  CZ  TYR A  79      27.481   2.320  65.126  1.00  0.00     0.115 A 
BRANCH 108 109
ATOM    109  OH  TYR A  79      28.181   1.360  65.840  1.00  0.00    -0.507 OA
ATOM    110  HH  TYR A  79      27.854   1.280  66.610  1.00  0.00     0.294 HD
ENDBRANCH 108 109
ENDBRANCH 102 103
ENDBRANCH 101 102
END_RES TYR A 79
BEGIN_RES THR A 222
ROOT
ATOM      1  CA  THR A 222      27.454   3.284  54.504  1.00  0.00     0.199 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  THR A 222      27.925   4.700  54.192  1.00  0.00     0.143 C 
ATOM      3 CG2  THR A 222      26.927   5.445  53.297  1.00  0.00     0.041 C 
BRANCH   2   4
ATOM      4 OG1  THR A 222      28.069   5.456  55.404  1.00  0.00    -0.389 OA
ATOM      5 HG1  THR A 222      28.805   5.278  55.768  1.00  0.00     0.211 HD
ENDBRANCH   2   4
ENDBRANCH   1   2
END_RES THR A 222
";

    /// The backbone N of each residue, which stays in the rigid file.
    const RIGID: &str = "\
ATOM    694  N   TYR A  79      24.786   6.752  64.683  1.00  0.00    -0.343 N 
ATOM   1954  N   THR A 222      28.499   2.529  55.204  1.00  0.00    -0.341 N 
";

    fn parse_flex_text(flex: &str) -> Result<FlexResidues, IoError> {
        let dir = tempfile::tempdir().unwrap();
        let (flex_path, rigid_path) = (
            dir.path().join("flex.pdbqt"),
            dir.path().join("rigid.pdbqt"),
        );
        std::fs::write(&flex_path, flex).unwrap();
        std::fs::write(&rigid_path, RIGID).unwrap();
        parse_flex_pdbqt(&flex_path, &parse_pdbqt(&rigid_path).unwrap())
    }

    /// Residues, torsions and the fixed atoms: each ROOT, and the first atom
    /// of each branch off it, which lies on that branch's axis.
    #[test]
    fn test_flex_residues() {
        let flex = parse_flex_text(FLEX).unwrap();
        let names: Vec<&str> = flex
            .molecule
            .atoms
            .iter()
            .map(|a| a.name.as_str())
            .collect();
        assert_eq!(flex.num_residues(), 2);
        assert_eq!(flex.residue_of.iter().filter(|&&r| r == 1).count(), 5);
        assert_eq!(flex.molecule.torsions.len(), 5);
        let fixed: Vec<&str> = (0..names.len())
            .filter(|&i| flex.fixed[i])
            .map(|i| names[i])
            .collect();
        assert_eq!(fixed, ["CA", "CB", "CA", "CB"]);
        // THR's chi1 (CA-CB, found by serial within its residue) turns the
        // atoms from CB on.
        let chi1 = flex.molecule.torsions.last().unwrap();
        let bond = &flex.molecule.bonds[chi1.bond_idx];
        assert_eq!((bond.atom1_idx, bond.atom2_idx), (10, 11));
        assert_eq!(chi1.moving_atoms, [11, 12, 13, 14]);
        // Each CA is bonded to its residue's rigid N, so it is polar.
        assert!(!flex.molecule.atoms[0].is_hydrophobic);
        assert!(flex.molecule.atoms[1].is_hydrophobic);
        assert!(flex.molecule.atoms[8].is_hbond_donor);
        assert!(flex
            .molecule
            .bonds
            .iter()
            .all(|b| b.atom1_idx < 15 && b.atom2_idx < 15));
    }

    /// Vina's flex pairs: heavy atoms more than three bonds apart that can
    /// move relative to each other. OH sits on the CZ-OH axis, so it is
    /// rigid relative to CB and CG, the ends of the ring branch's axis, but
    /// not to the CA.
    #[test]
    fn test_flex_internal_pairs() {
        let flex = parse_flex_text(FLEX).unwrap();
        let atoms = &flex.molecule.atoms;
        let label = |i: usize| format!("{}{}", atoms[i].name, atoms[i].residue_num);
        let within_tyr: Vec<String> = flex
            .internal_pairs()
            .into_iter()
            .filter(|&(i, j)| flex.residue_of[i] == 0 && flex.residue_of[j] == 0)
            .map(|(i, j)| format!("{}-{}", label(i), label(j)))
            .collect();
        assert_eq!(
            within_tyr,
            ["CA79-CE179", "CA79-CE279", "CA79-CZ79", "CA79-OH79"]
        );
        // Between residues, every heavy pair but the fixed ones.
        let between = flex
            .internal_pairs()
            .into_iter()
            .filter(|&(i, j)| flex.residue_of[i] != flex.residue_of[j])
            .count();
        let heavy = |r: usize| {
            (0..atoms.len())
                .filter(|&i| flex.residue_of[i] == r && !atoms[i].atom_type.is_hydrogen())
                .count()
        };
        assert_eq!(between, heavy(0) * heavy(1) - 4);
    }

    #[test]
    fn test_flex_errors() {
        let outside =
            "ATOM      1  CA  THR A 222      27.454   3.284  54.504  1.00  0.00     0.199 C \n";
        assert!(parse_flex_text(outside).is_err());
        let unclosed = FLEX.replace("END_RES THR A 222\n", "");
        assert!(parse_flex_text(&unclosed).is_err());
        let bad_branch = FLEX.replace("BRANCH   2   4", "BRANCH  99   4");
        assert!(parse_flex_text(&bad_branch).is_err());
    }
}
