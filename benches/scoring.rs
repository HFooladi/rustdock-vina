use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rustdock_vina::forcefield::vina::VinaForceField;
use rustdock_vina::molecule::Molecule;
use rustdock_vina::scoring::vina_score::VinaScore;

fn bench_vina_scoring(c: &mut Criterion) {
    // Create a small test molecule
    let mol = Molecule::new_test_molecule();
    let forcefield = VinaForceField::new();

    c.bench_function("vina_scoring", |b| {
        b.iter(|| {
            let score = VinaScore::new(&mol, &forcefield);
            black_box(score.calculate_score());
        })
    });
}

fn bench_ad4_scoring(c: &mut Criterion) {
    // Create a small test molecule
    let mol = Molecule::new_test_molecule();
    let forcefield = VinaForceField::new();

    c.bench_function("ad4_scoring", |b| {
        b.iter(|| {
            let score = VinaScore::new(&mol, &forcefield);
            black_box(score.calculate_ad4_score());
        })
    });
}

criterion_group!(scoring_benches, bench_vina_scoring, bench_ad4_scoring);
criterion_main!(scoring_benches);
