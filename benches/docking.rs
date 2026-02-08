use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nalgebra::Vector3;
use rustdock_vina::forcefield::vina::VinaForceField;
use rustdock_vina::io::parse_pdbqt;
use rustdock_vina::optimization::monte_carlo::{MonteCarlo, MonteCarloParams};
use rustdock_vina::optimization::Optimizer;
use std::path::PathBuf;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("test_data")
}

fn bench_complete_docking(c: &mut Criterion) {
    let receptor = parse_pdbqt(&test_data_dir().join("receptor.pdbqt")).unwrap();
    let ligand = parse_pdbqt(&test_data_dir().join("ligand.pdbqt")).unwrap();
    let forcefield = VinaForceField::new();

    let params = MonteCarloParams {
        steps_without_improvement: 50,
        use_local_optimization: false,
        ..MonteCarloParams::default()
    };
    let optimizer = MonteCarlo::with_params(params);

    let center = Vector3::new(5.0, 4.0, 2.0);
    let box_size = Vector3::new(10.0, 10.0, 10.0);

    c.bench_function("complete_docking", |b| {
        b.iter(|| {
            let results =
                optimizer.generate_poses(&ligand, &receptor, &forcefield, center, box_size, 2, 200);
            let _ = black_box(results);
        })
    });
}

criterion_group!(docking_benches, bench_complete_docking);
criterion_main!(docking_benches);
