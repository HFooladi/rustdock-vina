use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rustdock_vina::forcefield::vina::VinaForceField;
use rustdock_vina::molecule::Molecule;
use rustdock_vina::optimization::genetic::GeneticAlgorithm;
use rustdock_vina::optimization::monte_carlo::MonteCarloOptimizer;

fn bench_complete_docking_genetic(c: &mut Criterion) {
    let receptor = Molecule::new_test_receptor();
    let ligand = Molecule::new_test_ligand();
    let forcefield = VinaForceField::new();

    c.bench_function("complete_docking_genetic", |b| {
        b.iter(|| {
            let mut ga = GeneticAlgorithm::new(&receptor, &forcefield);
            let result = ga.dock(&ligand);
            black_box(result);
        })
    });
}

fn bench_complete_docking_monte_carlo(c: &mut Criterion) {
    let receptor = Molecule::new_test_receptor();
    let ligand = Molecule::new_test_ligand();
    let forcefield = VinaForceField::new();

    c.bench_function("complete_docking_monte_carlo", |b| {
        b.iter(|| {
            let mut mc = MonteCarloOptimizer::new(&receptor, &forcefield);
            let result = mc.dock(&ligand);
            black_box(result);
        })
    });
}

criterion_group!(
    docking_benches,
    bench_complete_docking_genetic,
    bench_complete_docking_monte_carlo
);
criterion_main!(docking_benches);
