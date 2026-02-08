use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rustdock_vina::forcefield::vina::VinaForceField;
use rustdock_vina::forcefield::ForceField;
use rustdock_vina::io::parse_pdbqt;
use std::path::PathBuf;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("test_data")
}

fn bench_vina_scoring(c: &mut Criterion) {
    let receptor = parse_pdbqt(&test_data_dir().join("receptor.pdbqt")).unwrap();
    let ligand = parse_pdbqt(&test_data_dir().join("ligand.pdbqt")).unwrap();
    let forcefield = VinaForceField::new();

    c.bench_function("vina_scoring", |b| {
        b.iter(|| {
            let mut total_energy = 0.0;
            for lig_atom in &ligand.atoms {
                for rec_atom in &receptor.atoms {
                    let distance = lig_atom.distance(rec_atom);
                    if distance < 8.0 && distance > 0.01 {
                        if let Ok(energy) =
                            forcefield.atom_pair_energy(lig_atom, rec_atom, distance)
                        {
                            total_energy += energy;
                        }
                    }
                }
            }
            black_box(total_energy);
        })
    });
}

criterion_group!(scoring_benches, bench_vina_scoring);
criterion_main!(scoring_benches);
