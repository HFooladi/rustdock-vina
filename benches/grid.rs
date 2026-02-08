use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nalgebra::Vector3;
use rustdock_vina::atom::AtomType;
use rustdock_vina::grid::Grid;

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

criterion_group!(grid_benches, bench_grid_creation, bench_grid_get_value);
criterion_main!(grid_benches);
