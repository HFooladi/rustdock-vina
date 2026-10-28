# Rustdock-Vina

A Rust implementation of [AutoDock Vina](https://vina.scripps.edu/) for molecular docking.

rustdock-vina reimplements AutoDock Vina 1.2.7's scoring and search. On the
CASF-2016 and PoseBusters benchmark sets it scores crystal poses the same as
Vina, docks about as accurately, and is faster, as the tables below show. Its
output files can be read and rescored by Vina.

## Features

- Vina and Vinardo scoring, matching Vina 1.2.7: the same atom typing, bond
  inference and torsional normalization
- Vina's Monte Carlo search: independent chains, BFGS local optimization with
  analytic gradients, and RMSD clustering of binding modes
- Precomputed per-atom-type affinity grids
- Batch screening: many ligands docked into one box in parallel, sharing one
  set of receptor grids
- PDBQT output in Vina's layout that keeps the input ligand's records, so
  Meeko and Vina read it back
- Reproducible runs with `--seed`
- Experimental AutoDock4 scoring (`--scoring ad4`; not validated)

## Parity with AutoDock Vina

Both programs docked the same prepared inputs: a Meeko ligand built from an
RDKit conformer, so the crystal pose never reaches the search, and an
Open Babel receptor. The search box is the crystal ligand's extent plus 10 Å.
Both used exhaustiveness 8, 9 modes and 16 CPUs, and every complex was docked
with five seeds (42, 1, 2, 3, 4). A top-1 success is a top pose within 2 Å
heavy-atom RMSD of the crystal ligand. Rates are mean ± standard deviation
over the seeds; a failed dock counts as a miss. The harness and the full
per-complex reports are in [`benchmark/parity/`](benchmark/parity/).

**Vina scoring** ([CASF-2016](benchmark/parity/PARITY_full_seeds.md),
[PoseBusters](benchmark/parity/PARITY_posebusters_seeds.md))

| | CASF-2016 (285): Vina | CASF-2016: rustdock | PoseBusters (428): Vina | PoseBusters: rustdock |
|---|---|---|---|---|
| Docked without error | 280 | 285 | 424 | 428 |
| Top-1 success (< 2 Å, %) | 53.3 ± 0.7 | 54.3 ± 0.2 | 54.7 ± 0.8 | 55.2 ± 0.3 |
| Any-mode success (< 2 Å, %) | 78.2 ± 0.7 | 79.2 ± 1.7 | 77.6 ± 1.4 | 79.0 ± 0.8 |
| Median top-1 RMSD (Å) | 1.65 | 1.61 | 1.77 | 1.77 |
| Top-1 PoseBusters-valid (%) | 90.7 | 92.8 | 93.4 | 93.7 |
| Median wall time (s) | 4.0 | 3.0 | 5.3 | 4.5 |
| Paired top-1 difference, rustdock − Vina (95% CI) | | −0.1 pp (−1.4 to +1.2) | | −0.2 pp (−1.6 to +1.1) |

**Vinardo scoring** ([CASF-2016](benchmark/parity/PARITY_full_vinardo_seeds.md),
[PoseBusters](benchmark/parity/PARITY_posebusters_vinardo_seeds.md))

| | CASF-2016 (285): Vina | CASF-2016: rustdock | PoseBusters (428): Vina | PoseBusters: rustdock |
|---|---|---|---|---|
| Docked without error | 280 | 285 | 424 | 428 |
| Top-1 success (< 2 Å, %) | 51.5 ± 1.1 | 54.0 ± 0.8 | 49.3 ± 0.6 | 51.6 ± 0.4 |
| Any-mode success (< 2 Å, %) | 72.2 ± 2.8 | 75.3 ± 1.3 | 71.2 ± 0.9 | 74.3 ± 0.8 |
| Median top-1 RMSD (Å) | 1.84 | 1.76 | 2.01 | 1.91 |
| Top-1 PoseBusters-valid (%) | 90.8 | 92.4 | 92.8 | 93.6 |
| Median wall time (s) | 4.3 | 3.3 | 5.7 | 4.6 |
| Paired top-1 difference, rustdock − Vina (95% CI) | | +2.2 pp (+0.6 to +3.9) | | +1.8 pp (+0.5 to +3.1) |

The paired difference compares each complex's top-1 success rate over the
five seeds, on the complexes both programs docked; the CI is a bootstrap over
complexes. Single-seed reports (seed 42) with every complex are also in
`benchmark/parity/` (`PARITY_full.md` and so on).

- **Scoring.** For all 704 crystal poses that Vina could read, rustdock's
  score is within 0.03 kcal/mol of Vina's, with either scoring function.
- **Search.** With Vina scoring, top-1 success is the same as Vina's to
  within about 1.5 percentage points either way. With Vinardo, rustdock finds
  the crystal pose about 2 points more often than Vina, a small difference
  that is unlikely to be chance (Wilcoxon p ≤ 0.02 on both sets). Its cause
  is not yet known.
- **Errors.** Vina rejects receptors containing As, Xe, Mo or B; rustdock
  docks them.
- **Output.** Vina reads rustdock's output files. It rescores the top pose
  at a median 0.002–0.003 kcal/mol from the affinity rustdock reports, and
  at most 0.020.
- **Screening.** Docking 100 CASF-2016 ligands into one pocket on 96 threads
  took 141 s with `--batch`, against 1024 s with Vina 1.2.7's `--batch` and
  210 s with 12 Vina processes of 8 threads. The best affinities correlate
  with Vina's at Pearson 0.98.

## Limitations

- **Flexible side chains are not implemented.** `--flex` is rejected with
  an error.
- **AD4 scoring is experimental.** It has not been validated against
  AutoDock4 and warns when selected.
- **PDBQT is the only input format.** There is no preparation step, so
  prepare inputs with Meeko, ADFR or MGLTools.
- **Config files** set only `center_x/y/z` and `size_x/y/z`; other keys are
  reported and ignored.

## Installation

### Prerequisites

- Rust 1.70 or later
- Cargo (Rust's package manager)

### Building from source

```bash
git clone https://github.com/HFooladi/rustdock-vina
cd rustdock-vina
cargo build --release
```

The compiled binary will be available at `target/release/rustdock`. It is
named `rustdock` rather than `vina` so that it can be installed next to
AutoDock Vina.

## Usage

### Docking

```bash
rustdock dock --receptor protein.pdbqt --ligand ligand.pdbqt --center 0,0,0 --size 20,20,20
```

All `dock` options:

| Flag | Description | Default |
|------|-------------|---------|
| `--receptor <FILE>` | PDBQT file containing the receptor (required) | |
| `--ligand <FILE>` | PDBQT file containing the ligand to dock | |
| `--flex <FILE>` | **Not implemented**; passing it is an error | |
| `-c, --config <FILE>` | Configuration file for docking parameters | |
| `--center <X,Y,Z>` | Center of the search box | |
| `--size <X,Y,Z>` | Size of the search box | |
| `--scoring <NAME>` | Scoring function: `vina`, `vinardo`, or experimental `ad4` | `vina` |
| `--exhaustiveness <N>` | Search exhaustiveness (higher = more accurate, slower) | `8` |
| `--num-modes <N>` | Number of binding modes to generate | `9` |
| `-o, --out <FILE>` | Output file for docking results | |
| `--batch <FILE/DIR>...` | Ligands to screen: PDBQT files, or directories of them | |
| `--dir <DIR>` | Output directory (required with `--batch`) | |
| `--cpu <N>` | Number of threads | all cores |
| `--energy-range <KCAL>` | Energy range for output poses (kcal/mol) | `3.0` |
| `--seed <N>` | Random seed; the same seed and inputs give the same poses | random |
| `--grid-spacing <A>` | Spacing of the precomputed affinity grids | `0.375` |
| `--no-grid` | Search on exact pairwise energies (slow) | |

### Screening a library

```bash
rustdock dock --receptor protein.pdbqt --batch ligands/ --dir results/ \
  --center 0,0,0 --size 20,20,20 --seed 42
```

`--batch` docks many ligands into one box. The receptor's grids are built
once and shared, and ligands are docked in parallel, so all cores stay busy
even though one dock runs only `exhaustiveness` chains. Each ligand's poses go
to `results/<name>_out.pdbqt`, the same as docking it alone with the same
seed. `results/scores.tsv` lists every ligand's best affinity, in input order.
A ligand that fails to parse or dock is listed there with its error, and the
run carries on.

### Scoring

Score an existing ligand pose against a receptor:

```bash
rustdock score --receptor protein.pdbqt --ligand ligand.pdbqt
```

All `score` options:

| Flag | Description | Default |
|------|-------------|---------|
| `--receptor <FILE>` | PDBQT file containing the receptor (required) | |
| `--ligand <FILE>` | PDBQT file containing the ligand pose (required) | |
| `--flex <FILE>` | **Not implemented**; passing it is an error | |
| `--scoring <NAME>` | Scoring function: `vina`, `vinardo`, or experimental `ad4` | `vina` |

### Configuration file

Instead of passing `--center` and `--size` on the command line, you can use a configuration file (`key = value` format, `#` for comments):

```ini
# Search box center
center_x = 15.190
center_y = 53.903
center_z = 16.917

# Search box size (Angstroms)
size_x = 20.0
size_y = 20.0
size_z = 20.0
```

```bash
rustdock dock --receptor protein.pdbqt --ligand ligand.pdbqt --config params.txt
```

### Example with test data

```bash
cargo run --bin rustdock -- dock \
  --receptor tests/test_data/receptor.pdbqt \
  --ligand tests/test_data/ligand.pdbqt \
  --center 0,0,0 --size 20,20,20

cargo run --bin rustdock -- score \
  --receptor tests/test_data/receptor.pdbqt \
  --ligand tests/test_data/ligand.pdbqt
```

## Project Structure

```
src/
  atom/           Atom types, van der Waals radii, H-bond properties
  molecule/       Molecule struct: atoms, bonds, torsions, geometry
  forcefield/     ForceField trait + Vina, Vinardo and AD4 implementations
  scoring/        High-level scoring combining forcefield terms
  optimization/   Optimizer trait + Vina's Monte Carlo search
  grid/           Precomputed affinity grids
  io/             PDBQT parser and writer
  math/           Vector/matrix utilities
  utils/          General helpers
  lib.rs          Library root
  main.rs         CLI binary
tests/            Integration tests and test data
benches/          Criterion benchmarks
benchmark/parity/ Parity harness and reports against AutoDock Vina
```

## Testing and Benchmarks

```bash
cargo test          # Run all tests
cargo bench         # Run performance benchmarks
```

## Documentation

Generate API documentation locally:

```bash
cargo doc --no-deps --document-private-items
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for development guidelines.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please read the [Contributing Guidelines](CONTRIBUTING.md) for details on submitting pull requests, reporting issues, and coding style.

## Acknowledgments

- [AutoDock Vina](https://vina.scripps.edu/) developers for the original implementation
- Trott, O., & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *Journal of Computational Chemistry*, 31(2), 455-461.
