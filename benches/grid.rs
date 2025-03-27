use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rustdock_vina::grid::Grid;
use rustdock_vina::molecule::Molecule;
use nalgebra::Point3;

fn bench_grid_creation(c: &mut Criterion) {
    c.bench_function("grid_creation", |b| {
        b.iter(|| {
            let grid = Grid::new(
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(20.0, 20.0, 20.0),
                0.375, // grid spacing
            );
            black_box(grid);
        })
    });
}

fn bench_grid_interpolation(c: &mut Criterion) {
    let grid = Grid::new(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(20.0, 20.0, 20.0),
        0.375,
    );
    
    c.bench_function("grid_interpolation", |b| {
        b.iter(|| {
            let point = Point3::new(10.0, 10.0, 10.0);
            black_box(grid.interpolate(&point));
        })
    });
}

fn bench_grid_population(c: &mut Criterion) {
    let mol = Molecule::new_test_molecule();
    let mut grid = Grid::new(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(20.0, 20.0, 20.0),
        0.375,
    );
    
    c.bench_function("grid_population", |b| {
        b.iter(|| {
            black_box(grid.populate_from_molecule(&mol));
        })
    });
}

criterion_group!(
    grid_benches,
    bench_grid_creation,
    bench_grid_interpolation,
    bench_grid_population
);
criterion_main!(grid_benches); 