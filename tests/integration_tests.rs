//! Integration tests for the rustdock-vina molecular docking library

use nalgebra::Vector3;
use rustdock_vina::forcefield::vina::VinaForceField;
use rustdock_vina::forcefield::ForceField;
use rustdock_vina::io::{parse_pdbqt, write_docking_results};
use rustdock_vina::optimization::monte_carlo::{MonteCarlo, MonteCarloParams};
use rustdock_vina::optimization::Optimizer;
use std::path::PathBuf;
use tempfile::tempdir;

/// Get the path to test data directory
fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("test_data")
}

#[test]
fn test_parse_receptor_pdbqt() {
    let receptor_path = test_data_dir().join("receptor.pdbqt");
    let receptor = parse_pdbqt(&receptor_path).expect("Failed to parse receptor PDBQT");

    assert!(!receptor.atoms.is_empty(), "Receptor should have atoms");
    assert!(
        receptor.atoms.len() >= 10,
        "Receptor should have at least 10 atoms"
    );
}

#[test]
fn test_parse_ligand_pdbqt() {
    let ligand_path = test_data_dir().join("ligand.pdbqt");
    let ligand = parse_pdbqt(&ligand_path).expect("Failed to parse ligand PDBQT");

    assert!(!ligand.atoms.is_empty(), "Ligand should have atoms");
    assert!(
        ligand.atoms.len() >= 4,
        "Ligand should have at least 4 atoms"
    );
    // Note: Torsion parsing depends on proper BRANCH record handling
    // For a simple test ligand, torsions may or may not be populated
}

#[test]
fn test_vina_scoring_basic() {
    let forcefield = VinaForceField::default();

    assert_eq!(forcefield.name(), "Vina");

    let receptor_path = test_data_dir().join("receptor.pdbqt");
    let ligand_path = test_data_dir().join("ligand.pdbqt");

    let receptor = parse_pdbqt(&receptor_path).expect("Failed to parse receptor");
    let ligand = parse_pdbqt(&ligand_path).expect("Failed to parse ligand");

    // Calculate energy between all atom pairs
    let mut total_energy = 0.0;
    let cutoff = 8.0;

    for lig_atom in &ligand.atoms {
        for rec_atom in &receptor.atoms {
            let distance = lig_atom.distance(rec_atom);
            if distance < cutoff && distance > 0.01 {
                if let Ok(energy) = forcefield.atom_pair_energy(lig_atom, rec_atom, distance) {
                    total_energy += energy;
                }
            }
        }
    }

    // Energy should be finite and reasonable
    assert!(total_energy.is_finite(), "Energy should be finite");
}

#[test]
fn test_monte_carlo_optimization() {
    let receptor_path = test_data_dir().join("receptor.pdbqt");
    let ligand_path = test_data_dir().join("ligand.pdbqt");

    let receptor = parse_pdbqt(&receptor_path).expect("Failed to parse receptor");
    let ligand = parse_pdbqt(&ligand_path).expect("Failed to parse ligand");

    let forcefield = VinaForceField::default();

    // Create optimizer with fewer iterations for faster testing
    let params = MonteCarloParams {
        steps_without_improvement: 100,
        use_local_optimization: false, // Disable for faster testing
        ..MonteCarloParams::default()
    };
    let optimizer = MonteCarlo::with_params(params);

    let center = Vector3::new(5.0, 4.0, 2.0);
    let box_size = Vector3::new(10.0, 10.0, 10.0);

    let result = optimizer
        .optimize(&ligand, &receptor, &forcefield, center, box_size, 500)
        .expect("Optimization should succeed");

    // Check that we got a valid result
    assert!(result.energy.is_finite(), "Energy should be finite");
    assert!(
        !result.molecule.atoms.is_empty(),
        "Result should have atoms"
    );
    assert!(
        !result.energy_components.is_empty(),
        "Should have energy components"
    );
}

#[test]
fn test_generate_multiple_poses() {
    let receptor_path = test_data_dir().join("receptor.pdbqt");
    let ligand_path = test_data_dir().join("ligand.pdbqt");

    let receptor = parse_pdbqt(&receptor_path).expect("Failed to parse receptor");
    let ligand = parse_pdbqt(&ligand_path).expect("Failed to parse ligand");

    let forcefield = VinaForceField::default();

    let params = MonteCarloParams {
        steps_without_improvement: 50,
        use_local_optimization: false,
        ..MonteCarloParams::default()
    };
    let optimizer = MonteCarlo::with_params(params);

    let center = Vector3::new(5.0, 4.0, 2.0);
    let box_size = Vector3::new(10.0, 10.0, 10.0);
    let num_poses = 3;

    let results = optimizer
        .generate_poses(
            &ligand,
            &receptor,
            &forcefield,
            center,
            box_size,
            num_poses,
            200,
        )
        .expect("Should generate poses");

    assert_eq!(
        results.len(),
        num_poses,
        "Should generate requested number of poses"
    );

    // Results should be sorted by energy
    for i in 1..results.len() {
        assert!(
            results[i].energy >= results[i - 1].energy,
            "Results should be sorted by energy"
        );
    }

    // All but the first should have RMSD
    for result in results.iter().skip(1) {
        assert!(result.rmsd.is_some(), "Non-first poses should have RMSD");
    }
}

#[test]
fn test_write_docking_results() {
    let receptor_path = test_data_dir().join("receptor.pdbqt");
    let ligand_path = test_data_dir().join("ligand.pdbqt");

    let receptor = parse_pdbqt(&receptor_path).expect("Failed to parse receptor");
    let ligand = parse_pdbqt(&ligand_path).expect("Failed to parse ligand");

    let forcefield = VinaForceField::default();

    let params = MonteCarloParams {
        steps_without_improvement: 50,
        use_local_optimization: false,
        ..MonteCarloParams::default()
    };
    let optimizer = MonteCarlo::with_params(params);

    let center = Vector3::new(5.0, 4.0, 2.0);
    let box_size = Vector3::new(10.0, 10.0, 10.0);

    let results = optimizer
        .generate_poses(&ligand, &receptor, &forcefield, center, box_size, 2, 100)
        .expect("Should generate poses");

    // Write results to a temporary file
    let temp_dir = tempdir().expect("Failed to create temp dir");
    let output_path = temp_dir.path().join("output.pdbqt");

    write_docking_results(&results, &output_path).expect("Should write results");

    // Verify the file was created and has content
    assert!(output_path.exists(), "Output file should exist");

    let content = std::fs::read_to_string(&output_path).expect("Should read output file");
    assert!(!content.is_empty(), "Output file should not be empty");
    assert!(
        content.contains("MODEL"),
        "Output should contain MODEL records"
    );
    assert!(
        content.contains("ATOM"),
        "Output should contain ATOM records"
    );
}

#[test]
fn test_forcefield_individual_terms() {
    let forcefield = VinaForceField::default();

    let receptor_path = test_data_dir().join("receptor.pdbqt");
    let ligand_path = test_data_dir().join("ligand.pdbqt");

    let receptor = parse_pdbqt(&receptor_path).expect("Failed to parse receptor");
    let ligand = parse_pdbqt(&ligand_path).expect("Failed to parse ligand");

    // Test individual energy terms
    if let (Some(lig_atom), Some(rec_atom)) = (ligand.atoms.first(), receptor.atoms.first()) {
        let distance = lig_atom.distance(rec_atom);

        if distance < 8.0 && distance > 0.01 {
            let vdw = forcefield.vdw_energy(lig_atom, rec_atom, distance);
            let elec = forcefield.electrostatic_energy(lig_atom, rec_atom, distance);
            let hbond = forcefield.hbond_energy(lig_atom, rec_atom, distance);
            let hydrophobic = forcefield.hydrophobic_energy(lig_atom, rec_atom, distance);
            let desolv = forcefield.desolvation_energy(lig_atom, rec_atom, distance);

            // All terms should be computable
            assert!(vdw.is_ok(), "vdW energy should be computable");
            assert!(elec.is_ok(), "Electrostatic energy should be computable");
            assert!(hbond.is_ok(), "H-bond energy should be computable");
            assert!(
                hydrophobic.is_ok(),
                "Hydrophobic energy should be computable"
            );
            assert!(desolv.is_ok(), "Desolvation energy should be computable");

            // All energies should be finite
            assert!(vdw.unwrap().is_finite());
            assert!(elec.unwrap().is_finite());
            assert!(hbond.unwrap().is_finite());
            assert!(hydrophobic.unwrap().is_finite());
            assert!(desolv.unwrap().is_finite());
        }
    }
}

#[test]
fn test_monte_carlo_seed_is_deterministic() {
    let receptor = parse_pdbqt(test_data_dir().join("receptor.pdbqt")).unwrap();
    let ligand = parse_pdbqt(test_data_dir().join("ligand.pdbqt")).unwrap();
    let forcefield = VinaForceField::default();

    let params = MonteCarloParams {
        steps_without_improvement: 50,
        use_local_optimization: false,
        seed: Some(42),
        ..MonteCarloParams::default()
    };
    let optimizer = MonteCarlo::with_params(params);

    let center = Vector3::new(5.0, 4.0, 2.0);
    let box_size = Vector3::new(10.0, 10.0, 10.0);

    let a = optimizer
        .optimize(&ligand, &receptor, &forcefield, center, box_size, 200)
        .unwrap();
    let b = optimizer
        .optimize(&ligand, &receptor, &forcefield, center, box_size, 200)
        .unwrap();

    assert_eq!(a.energy, b.energy, "same seed must yield same energy");
    assert_eq!(
        a.molecule.atoms.len(),
        b.molecule.atoms.len(),
        "same seed must yield same atom count"
    );
    for (atom_a, atom_b) in a.molecule.atoms.iter().zip(&b.molecule.atoms) {
        assert_eq!(
            atom_a.coordinates, atom_b.coordinates,
            "same seed must yield identical pose coordinates"
        );
    }
}

#[test]
fn test_generate_poses_seed_is_deterministic() {
    let receptor = parse_pdbqt(test_data_dir().join("receptor.pdbqt")).unwrap();
    let ligand = parse_pdbqt(test_data_dir().join("ligand.pdbqt")).unwrap();
    let forcefield = VinaForceField::default();

    let params = MonteCarloParams {
        steps_without_improvement: 50,
        use_local_optimization: false,
        seed: Some(7),
        ..MonteCarloParams::default()
    };
    let optimizer = MonteCarlo::with_params(params);

    let center = Vector3::new(5.0, 4.0, 2.0);
    let box_size = Vector3::new(10.0, 10.0, 10.0);

    let first = optimizer
        .generate_poses(&ligand, &receptor, &forcefield, center, box_size, 3, 100)
        .unwrap();
    let second = optimizer
        .generate_poses(&ligand, &receptor, &forcefield, center, box_size, 3, 100)
        .unwrap();

    assert_eq!(first.len(), second.len());
    for (a, b) in first.iter().zip(&second) {
        assert_eq!(
            a.energy, b.energy,
            "same master seed must give same energies"
        );
    }
}

#[test]
fn test_molecule_center_and_bounding_box() {
    let ligand_path = test_data_dir().join("ligand.pdbqt");
    let ligand = parse_pdbqt(&ligand_path).expect("Failed to parse ligand");

    let center = ligand.center().expect("Should compute center");
    let (min_coord, max_coord) = ligand.bounding_box().expect("Should compute bounding box");

    // Center should be within bounding box
    assert!(center.x >= min_coord.x && center.x <= max_coord.x);
    assert!(center.y >= min_coord.y && center.y <= max_coord.y);
    assert!(center.z >= min_coord.z && center.z <= max_coord.z);
}

/// Path to the 1IEP benchmark fixtures (imatinib bound to Abl kinase).
fn benchmark_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("benchmark_data")
}

/// Scoring the crystal pose of 1IEP must land close to AutoDock Vina.
///
/// Real Vina reports an intermolecular energy of -16.026 for this complex,
/// which becomes -11.865 kcal/mol after dividing by `1 + 0.05846 * 6`. This is
/// the test that pins the atom typing: treating every carbon as hydrophobic
/// and inferring H-bond donors from `NA`/`OA` labels put this near -16.8,
/// roughly 5 kcal/mol too favourable, and nothing in the suite noticed.
#[test]
fn test_1iep_crystal_pose_scores_close_to_vina() {
    let receptor =
        parse_pdbqt(benchmark_data_dir().join("receptor.pdbqt")).expect("parse receptor");
    let ligand = parse_pdbqt(benchmark_data_dir().join("ligand.pdbqt")).expect("parse ligand");
    let forcefield = VinaForceField::default();

    let mut inter = 0.0;
    for lig_atom in &ligand.atoms {
        for rec_atom in &receptor.atoms {
            let d = lig_atom.distance(rec_atom);
            if d > 0.01 {
                inter += forcefield.atom_pair_energy(lig_atom, rec_atom, d).unwrap();
            }
        }
    }
    let affinity = forcefield.calculate_affinity(inter, ligand.torsions.len());

    const VINA_REFERENCE: f64 = -11.865;
    const TOLERANCE: f64 = 1.5;
    assert!(
        (affinity - VINA_REFERENCE).abs() < TOLERANCE,
        "1IEP crystal-pose affinity {affinity:.3} kcal/mol is more than \
         {TOLERANCE} from Vina's {VINA_REFERENCE}"
    );
}

/// The ligand's full connectivity must be inferred, not just the rotatable
/// bonds from the `BRANCH` records.
///
/// Imatinib has 37 heavy atoms and 6 rotatable bonds. Parsing used to stop at
/// those 6 bonds, which left every bond-graph-derived property wrong.
#[test]
fn test_ligand_connectivity_is_inferred() {
    let ligand = parse_pdbqt(benchmark_data_dir().join("ligand.pdbqt")).expect("parse ligand");

    assert_eq!(ligand.torsions.len(), 6, "imatinib has 6 rotatable bonds");
    assert!(
        ligand.bonds.len() >= ligand.atoms.len(),
        "expected full connectivity (>= {} bonds), got {}",
        ligand.atoms.len(),
        ligand.bonds.len()
    );
    assert!(
        !ligand.intramolecular_pairs().is_empty(),
        "a 37-atom ligand must have intramolecular scoring pairs"
    );
}

/// A carbon bonded to a heteroatom is `C_P` and must not be hydrophobic.
#[test]
fn test_polar_carbons_are_typed_non_hydrophobic() {
    let receptor =
        parse_pdbqt(benchmark_data_dir().join("receptor.pdbqt")).expect("parse receptor");

    let carbons = receptor
        .atoms
        .iter()
        .filter(|a| a.atom_type == rustdock_vina::AtomType::Carbon)
        .count();
    let hydrophobic_carbons = receptor
        .atoms
        .iter()
        .filter(|a| a.atom_type == rustdock_vina::AtomType::Carbon && a.is_hydrophobic)
        .count();

    assert!(carbons > 0);
    assert!(
        hydrophobic_carbons < carbons,
        "some carbons must be C_P; got {hydrophobic_carbons} of {carbons} hydrophobic"
    );
}

/// This receptor carries no polar hydrogens, so under Vina's typing it has no
/// hydrogen-bond donors at all.
#[test]
fn test_no_donors_without_polar_hydrogens() {
    let receptor =
        parse_pdbqt(benchmark_data_dir().join("receptor.pdbqt")).expect("parse receptor");

    assert!(
        !receptor.atoms.iter().any(|a| a.atom_type.is_hydrogen()),
        "fixture is expected to have no hydrogens"
    );
    assert!(
        !receptor.atoms.iter().any(|a| a.is_hbond_donor),
        "no atom can donate a hydrogen bond without a bonded polar hydrogen"
    );
}
