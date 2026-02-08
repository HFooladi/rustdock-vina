# CLAUDE.md

## Quick Reference

```bash
# Build & test (run these before committing)
cargo fmt                              # Format code
cargo clippy -- -D warnings            # Lint (must pass with zero warnings)
cargo test                             # Run all tests (unit + integration)
cargo build --release                  # Optimized build

# Run the CLI
cargo run --bin vina -- dock --receptor tests/test_data/receptor.pdbqt --ligand tests/test_data/ligand.pdbqt --center 0,0,0 --size 20,20,20
cargo run --bin vina -- score --receptor tests/test_data/receptor.pdbqt --ligand tests/test_data/ligand.pdbqt

# Other
cargo bench                            # Performance benchmarks (benches/{docking,grid,optimization,scoring}.rs)
cargo doc --no-deps --document-private-items  # Generate docs

# GitHub Actions
gh run list --limit 5                  # Check recent CI runs
gh run view <run-id> --log-failed      # Inspect CI failures
```

## CI Pipeline (.github/workflows/ci.yml)

CI runs on every push/PR to `main`. It has three jobs:

| Job | Steps | Notes |
|---|---|---|
| **Test** | `cargo fmt --all -- --check` → `cargo clippy -- -D warnings` → `cargo test --verbose` → `cargo bench --verbose` | Formatting is checked first; fix with `cargo fmt` before pushing |
| **Security** | `cargo install cargo-audit` → `cargo audit` | Checks for known vulnerabilities in dependencies |
| **Documentation** | `cargo doc --no-deps --document-private-items` | Ensures docs compile cleanly |

All three jobs are passing as of 2026-02-08.

**Before any commit, always run:** `cargo fmt && cargo clippy -- -D warnings && cargo test`

## Project Structure

Rust implementation of the AutoDock Vina molecular docking algorithm.

```
src/
├── lib.rs              # Library root — re-exports Atom, Molecule, VinaScore
├── main.rs             # CLI binary ("vina") — subcommands: dock, score
├── atom/               # Atom types (AtomType enum), van der Waals radii, H-bond properties
├── molecule/           # Molecule struct: atoms, bonds, torsions, geometric operations
├── forcefield/         # ForceField trait + implementations:
│   ├── vina.rs         #   VinaForceField (gauss, repulsion, hydrophobic, H-bond terms)
│   └── ad4.rs          #   AD4ForceField (AutoDock4 scoring)
├── scoring/            # VinaScore — high-level scoring combining forcefield terms
├── optimization/       # Optimizer trait + MonteCarlo implementation
│   └── monte_carlo.rs  #   Pose generation and BFGS-like optimization
├── grid/               # Spatial grid for efficient energy lookup
├── io/                 # PDBQT parser (parse_pdbqt) and writer (write_docking_results)
├── math/               # Vector/matrix utilities, coordinate transforms
└── utils/              # General helpers
tests/
├── integration_tests.rs  # 8 integration tests covering parsing, scoring, optimization, output
└── test_data/            # receptor.pdbqt, ligand.pdbqt
benches/                  # Criterion benchmarks: docking, grid, optimization, scoring
```

### Key Types

- `Molecule` — central struct for ligands/receptors (atoms, bonds, torsion tree)
- `Atom` — position, type (`AtomType` enum), charge, residue info
- `ForceField` trait — `calculate_energy(&self, ligand, receptor) -> f64`
- `Optimizer` trait — `generate_poses(ligand, receptor, ...) -> Vec<DockingResult>`
- `DockingResult` — energy, conformation, per-term breakdown

### Execution Flow

1. Parse PDBQT files → `Molecule` structs (`src/io/`)
2. Select force field (Vina or AD4) (`src/forcefield/`)
3. Monte Carlo optimization generates + optimizes poses (`src/optimization/`)
4. Score each pose with force field terms (`src/scoring/`)
5. Write results to PDBQT (`src/io/`)

## Conventions

- Library crate name: `rustdock_vina` (underscore). Binary name: `vina`.
- Unit tests live inside modules (`#[cfg(test)]`). Integration tests in `tests/`.
- All code changes must pass `cargo fmt`, `cargo clippy -- -D warnings`, and `cargo test` before merging.
- The release profile uses LTO, single codegen unit, and `opt-level = 3`.
