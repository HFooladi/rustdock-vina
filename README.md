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
- Flexible receptor side chains (`--flex`, Vina's flex residue PDBQT),
  searched and scored as Vina does
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
| Top-1 success (< 2 Å, %) | 53.3 ± 0.7 | 53.3 ± 0.8 | 54.7 ± 0.8 | 55.7 ± 0.7 |
| Any-mode success (< 2 Å, %) | 78.2 ± 0.7 | 77.6 ± 1.6 | 77.6 ± 1.4 | 78.6 ± 0.5 |
| Median top-1 RMSD (Å) | 1.65 | 1.74 | 1.77 | 1.76 |
| Top-1 PoseBusters-valid (%) | 90.7 | 92.6 | 93.4 | 94.2 |
| Median wall time (s) | 4.0 | 3.1 | 5.3 | 4.5 |
| Paired top-1 difference, rustdock − Vina (95% CI) | | −1.1 pp (−2.6 to +0.3) | | +0.3 pp (−1.0 to +1.6) |

**Vinardo scoring** ([CASF-2016](benchmark/parity/PARITY_full_vinardo_seeds.md),
[PoseBusters](benchmark/parity/PARITY_posebusters_vinardo_seeds.md))

| | CASF-2016 (285): Vina | CASF-2016: rustdock | PoseBusters (428): Vina | PoseBusters: rustdock |
|---|---|---|---|---|
| Docked without error | 280 | 285 | 424 | 428 |
| Top-1 success (< 2 Å, %) | 51.5 ± 1.1 | 52.9 ± 1.2 | 49.3 ± 0.6 | 49.5 ± 0.8 |
| Any-mode success (< 2 Å, %) | 72.2 ± 2.8 | 73.5 ± 0.9 | 71.2 ± 0.9 | 70.6 ± 0.6 |
| Median top-1 RMSD (Å) | 1.84 | 1.82 | 2.01 | 2.04 |
| Top-1 PoseBusters-valid (%) | 90.8 | 92.0 | 92.8 | 93.1 |
| Median wall time (s) | 4.3 | 3.6 | 5.7 | 4.8 |
| Paired top-1 difference, rustdock − Vina (95% CI) | | +0.8 pp (−0.5 to +2.1) | | −0.3 pp (−1.5 to +0.8) |

**Flexible side chains, Vina scoring** ([CASF-2016](benchmark/parity/PARITY_full_flex_seeds.md))

Up to three residues whose side chains are within 4 Å of the crystal
ligand were made flexible with Meeko's `mk_prepare_receptor.py --flexres`,
and the box was widened to cover them. Meeko could not prepare 51 of the 285
complexes, mostly ones with ions or metals it has no covalent radius for
([list](benchmark/parity/PARITY_full_flex_excluded.tsv)); they are left out
for both programs.

| | CASF-2016 (234): Vina | CASF-2016: rustdock |
|---|---|---|
| Top-1 success (< 2 Å, %) | 44.0 ± 1.2 | 44.4 ± 0.9 |
| Any-mode success (< 2 Å, %) | 72.5 ± 1.3 | 72.9 ± 1.5 |
| Median top-1 RMSD (Å) | 2.67 | 2.63 |
| Top-1 PoseBusters-valid (%) | 70.0 | 69.5 |
| Median wall time (s) | 15.6 | 10.6 |
| Paired top-1 difference, rustdock − Vina (95% CI) | | +0.4 pp (−1.1 to +2.1) |

PoseBusters checks the ligand against the protein with its crystal side
chains, so moved side chains it overlaps count against both programs.

**Virtual screening, Vina scoring** ([DUDE-Z](benchmark/screening/SCREENING.md))

Ten DUDE-Z targets (AA2AR, ABL1, ACES, AMPC, ANDR, EGFR, HIVPR, HS90A, PARP1,
THRB) with all their actives and decoys: 665 actives and 40,443 decoys,
docked by both programs from the same inputs with seeds 42 and 1, one thread
per dock. A compound scores its best top-pose affinity over its protomers.
The harness and protocol are in [`benchmark/screening/`](benchmark/screening/).

| Mean over 10 targets (seed 42) | Vina | rustdock | rustdock − Vina (95% CI) | Vina − Vina, two seeds |
|---|---|---|---|---|
| ROC AUC | 0.678 | 0.678 | +0.000 (−0.003 to +0.003) | +0.002 |
| BEDROC (α = 20) | 0.192 | 0.189 | −0.003 (−0.008 to +0.001) | −0.001 |
| EF1% | 7.1 | 6.7 | −0.37 (−1.21 to +0.69) | −0.22 |
| EF5% | 3.76 | 3.58 | −0.19 (−0.38 to +0.09) | −0.09 |
| Mean seconds per dock (one thread) | 26.7 | 23.8 | | |

The paired difference compares each complex's top-1 success rate over the
five seeds, on the complexes both programs docked; the CI is a bootstrap over
complexes. Single-seed reports (seed 42) with every complex are also in
`benchmark/parity/` (`PARITY_full.md` and so on).

- **Scoring.** For all 704 crystal poses that Vina could read, rustdock's
  score is within 0.03 kcal/mol of Vina's, with either scoring function.
- **Search.** Top-1 success is the same as Vina's to within about 1.5
  percentage points on every set, with both scoring functions and with
  flexible side chains; every 95% CI includes zero. On CASF-2016 with Vina
  scoring the rates are equal (53.3%), but on the 280 complexes Vina can
  dock rustdock is 1.1 points lower, a gap at the edge of chance (Wilcoxon
  p = 0.02, CI −2.6 to +0.3) that PoseBusters does not show. The docks only
  one program solves split both ways: on CASF-2016, 39 for Vina and 23 for
  rustdock; on PoseBusters, 70 and 76; together, 109 and 99. In most of them
  the loser's top pose scores worse than the winner's, so its search missed
  the crystal basin the other found, rather than scoring pulling it
  elsewhere. Top affinities have a median difference within 0.01 kcal/mol
  of Vina's.
- **The Vinardo edge was the search cap.** An earlier run had rustdock about
  2 points ahead of Vina with Vinardo. During Monte Carlo, rustdock capped
  each atom's energy against the receptor at 1.5 kcal/mol, the default in
  Vina's `monte_carlo.h`, where Vina 1.2.7's docking sets it to 10. With the
  cap at 10 the edge is gone.
- **Errors.** Vina rejects receptors containing As, Xe, Mo or B; rustdock
  docks them.
- **Output.** Vina reads rustdock's output files. It rescores the top pose
  at a median 0.002–0.004 kcal/mol from the affinity rustdock reports, and
  at most 0.019.
- **Virtual screening.** rustdock ranks actives above decoys as Vina does:
  every difference is inside its CI around zero and the size of Vina's own
  seed-to-seed difference. Per compound, its scores correlate with Vina's at
  Spearman 0.994, as Vina's do with themselves across seeds. rustdock's top
  affinity is a little lower more often than Vina's (19.9% vs 11.1% of
  compounds by more than 0.05 kcal/mol), and Vina's scoring agrees. That
  comes mostly from rustdock refining poses on exact pair energies, where
  Vina interpolates lookup tables
  ([analysis](benchmark/screening/README.md#rustdock-finds-slightly-lower-minima)).
- **Batch speed.** Docking 100 CASF-2016 ligands into one pocket on 96 threads
  took 141 s with `--batch`, against 1024 s with Vina 1.2.7's `--batch` and
  210 s with 12 Vina processes of 8 threads. The best affinities correlate
  with Vina's at Pearson 0.98.

## Limitations

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
| `--flex <FILE>` | Flexible side chains, as for Vina's `--flex`; `--receptor` is then the rigid part | |
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

### Flexible side chains

```bash
rustdock dock --receptor protein_rigid.pdbqt --flex protein_flex.pdbqt \
  --ligand ligand.pdbqt --center 0,0,0 --size 20,20,20 -o out.pdbqt
```

The side chains come in Vina's two files: the flexible residues as
`BEGIN_RES`…`END_RES` blocks, and the rigid receptor without them, as
Meeko's `mk_prepare_receptor.py --flexres` (or ADT's
`prepare_flexreceptor4.py`) writes them. As in Vina, the search turns each
residue's torsions along with the ligand's, and each mode in the output
file has the moved residues after the ligand. The side chains' energy with
the ligand counts toward the affinity; their energy with the rest of the
receptor and with each other is internal and only ranks the modes. The box
must cover the side chains as well as the ligand.

Docking CASF-2016 with flexible side chains, rustdock matches Vina's
success rate and takes about two thirds of its time (see
[Parity](#parity-with-autodock-vina)).

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
| `--flex <FILE>` | Flexible side chains, as for Vina's `--flex`; `--receptor` is then the rigid part | |
| `--scoring <NAME>` | Scoring function: `vina`, `vinardo`, or experimental `ad4` | `vina` |

With `--flex`, the breakdown adds the ligand's energy with the side chains,
which is part of the intermolecular energy and the affinity, and the side
chains' energies with the rigid receptor and with each other. As in Vina,
those two are internal energies and do not enter the affinity. Each
residue's ROOT atoms and the first atom of each branch off the ROOT (the CA
and CB in Meeko's and ADT's files) stay fixed and are not scored against the
receptor or the ligand.

Vina's own `--score_only --flex` affinity differs from rustdock's. Vina
subtracts an "unbound" energy that takes the side chains' energy with the
receptor from its grid, while the internal energy it adds back takes it
from exact pairs, so the two do not cancel. Passing Vina
`--unbound_energy` equal to its printed internal energy (2) gives
rustdock's affinity.

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
