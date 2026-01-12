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
