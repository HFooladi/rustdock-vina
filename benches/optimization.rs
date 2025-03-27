use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rustdock_vina::forcefield::vina::VinaForceField;
use rustdock_vina::molecule::Molecule;
use rustdock_vina::optimization::genetic::GeneticAlgorithm;
use rustdock_vina::optimization::monte_carlo::MonteCarloOptimizer;

fn bench_genetic_algorithm(c: &mut Criterion) {
    let mol = Molecule::new_test_molecule();
    let forcefield = VinaForceField::new();

    c.bench_function("genetic_algorithm", |b| {
        b.iter(|| {
            let mut ga = GeneticAlgorithm::new(&mol, &forcefield);
            black_box(ga.optimize(100)); // Run for 100 generations
        })
    });
}

fn bench_monte_carlo(c: &mut Criterion) {
    let mol = Molecule::new_test_molecule();
    let forcefield = VinaForceField::new();

    c.bench_function("monte_carlo", |b| {
        b.iter(|| {
            let mut mc = MonteCarloOptimizer::new(&mol, &forcefield);
            black_box(mc.optimize(1000)); // Run for 1000 iterations
        })
    });
}

criterion_group!(
    optimization_benches,
    bench_genetic_algorithm,
    bench_monte_carlo
);
criterion_main!(optimization_benches);
