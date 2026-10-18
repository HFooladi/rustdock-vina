use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nalgebra::Vector3;
use rustdock_vina::atom::AtomType;
use rustdock_vina::forcefield::vina::VinaForceField;
use rustdock_vina::grid::{Grid, GridSet};
use rustdock_vina::io::parse_pdbqt;
use std::path::PathBuf;

fn bench_grid_creation(c: &mut Criterion) {
    c.bench_function("grid_creation", |b| {
        b.iter(|| {
            let grid = Grid::new(
                Vector3::new(0.0, 0.0, 0.0),
                0.375,
                Vector3::new(54, 54, 54),
                AtomType::Carbon,
            );
            let _ = black_box(grid);
        })
    });
}

fn bench_grid_get_value(c: &mut Criterion) {
    let grid = Grid::new(
        Vector3::new(0.0, 0.0, 0.0),
        0.375,
        Vector3::new(54, 54, 54),
        AtomType::Carbon,
    )
    .expect("Failed to create grid");

    c.bench_function("grid_get_value", |b| {
        b.iter(|| {
            let point = Vector3::new(10.0, 10.0, 10.0);
            let _ = black_box(grid.get_value(&point));
        })
    });
}

/// Affinity grids for imatinib in Abl kinase (1IEP): 2,200 receptor heavy
/// atoms, a 20 A box plus margin at 0.375 A, as a real dock builds them.
fn bench_gridset_build(c: &mut Criterion) {
    let data = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("benchmark_data");
    let receptor = parse_pdbqt(data.join("receptor.pdbqt")).unwrap();
    let ligand = parse_pdbqt(data.join("ligand.pdbqt")).unwrap();
    let forcefield = VinaForceField::new();
    let center = ligand.center().unwrap();
    let box_size = Vector3::repeat(20.0);

    let mut group = c.benchmark_group("gridset");
    group.sample_size(10);
    group.bench_function("build_1iep", |b| {
        b.iter(|| {
            let grids = GridSet::build(
                &ligand,
                &receptor,
                &forcefield,
                &center,
                &box_size,
                4.0,
                0.375,
            );
            let _ = black_box(grids);
        })
    });
    group.finish();
}

criterion_group!(
    grid_benches,
    bench_grid_creation,
    bench_grid_get_value,
    bench_gridset_build
);
criterion_main!(grid_benches);
