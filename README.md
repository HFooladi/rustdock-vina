# Rustdock-Vina

A Rust implementation of [AutoDock Vina](https://vina.scripps.edu/) for molecular docking.

> **Status: work in progress — not yet validated for production use.**
> The scoring function now agrees with AutoDock Vina to about 1 kcal/mol on our
> one benchmark complex, but the conformational search is much weaker than
> Vina's and finds the correct pose only part of the time. Please do not treat
> results from this tool as a substitute for AutoDock Vina yet. See
> [Status and limitations](#status-and-limitations) for measured numbers.

## Features

- Vina scoring function (gauss, repulsion, hydrophobic, H-bond terms) with
  bond-graph-derived `xs` atom typing
- Monte Carlo search with local minimization, parallelized across poses via
  [rayon](https://github.com/rayon-rs/rayon)
- RMSD-clustered binding modes
- Configuration file support for the search box
- Batch ligand docking with directory output
- PDBQT file parsing and writing
- Memory-safe implementation in Rust

## Status and limitations

Measured on 1IEP (imatinib bound to Abl kinase), the one complex in
`benchmark_data/`, against real AutoDock Vina run on the same input files:

| | AutoDock Vina | rustdock-vina |
|---|---|---|
| Affinity of the crystal pose | −11.87 kcal/mol | **−10.93 kcal/mol** |
| Redocking success (best pose within 2 Å) | ~always | **~40% of runs** |
| Wall time (8 exhaustiveness, 24 cores) | ~1 s | ~14 s |

What this means in practice:

- **Scoring is close but not exact.** Roughly 1 kcal/mol high on this complex.
  A single complex is not a validation set.
- **The search is the weak part.** Runs that fail converge on a secondary site
  about 12 Å away at roughly −9 kcal/mol. Vina finds that site too and ranks it
  below the true one; the difference is that Vina reliably also finds the true
  one. Raising `--exhaustiveness` helps.
- **Performance is far off.** Every energy evaluation is a full O(N_ligand ×
  N_receptor) loop. There is no precomputed affinity grid, which is the single
  biggest reason for the gap.

Known gaps, all of which are unimplemented rather than partially working:

- **Vinardo scoring and flexible side chains are not implemented.** Earlier
  versions of this README advertised both; `--flex` is accepted and ignored.
- **The AD4 scoring function is experimental** and is not validated against
  AutoDock4 — it warns when selected. Its parameter table is a stub.
- **Output PDBQT is not round-trippable through AutoDock Vina**: no
  `ROOT`/`BRANCH`/`TORSDOF` records are written, and atom-name columns are
  offset by one, which also affects some structure viewers.
- **PDBQT is the only supported format.** There is no receptor preparation
  step, so inputs must come from MGLTools/ADFR/Meeko.
- Only `center_x/y/z` and `size_x/y/z` are read from a config file; any other
  key is reported and ignored.

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

The compiled binary will be available at `target/release/vina`.

## Usage

### Docking

```bash
vina dock --receptor protein.pdbqt --ligand ligand.pdbqt --center 0,0,0 --size 20,20,20
```

All `dock` options:

| Flag | Description | Default |
|------|-------------|---------|
| `--receptor <FILE>` | PDBQT file containing the receptor (required) | |
| `--ligand <FILE>` | PDBQT file containing the ligand to dock | |
| `--flex <FILE>` | Accepted but **not implemented** — currently ignored | |
| `-c, --config <FILE>` | Configuration file for docking parameters | |
| `--center <X,Y,Z>` | Center of the search box | |
| `--size <X,Y,Z>` | Size of the search box | |
| `--scoring <NAME>` | Scoring function: `vina`, or experimental `ad4` | `vina` |
| `--exhaustiveness <N>` | Search exhaustiveness (higher = more accurate, slower) | `8` |
| `--num-modes <N>` | Number of binding modes to generate | `9` |
| `-o, --out <FILE>` | Output file for docking results | |
| `--dir <DIR>` | Output directory for batch mode results | |
| `--energy-range <KCAL>` | Energy range for output poses (kcal/mol) | `3.0` |

### Scoring

Score an existing ligand pose against a receptor:

```bash
vina score --receptor protein.pdbqt --ligand ligand.pdbqt
```

All `score` options:

| Flag | Description | Default |
|------|-------------|---------|
| `--receptor <FILE>` | PDBQT file containing the receptor (required) | |
| `--ligand <FILE>` | PDBQT file containing the ligand pose (required) | |
| `--flex <FILE>` | Accepted but **not implemented** — currently ignored | |
| `--scoring <NAME>` | Scoring function: `vina`, or experimental `ad4` | `vina` |

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
vina dock --receptor protein.pdbqt --ligand ligand.pdbqt --config params.txt
```

### Example with test data

```bash
cargo run --bin vina -- dock \
  --receptor tests/test_data/receptor.pdbqt \
  --ligand tests/test_data/ligand.pdbqt \
  --center 0,0,0 --size 20,20,20

cargo run --bin vina -- score \
  --receptor tests/test_data/receptor.pdbqt \
  --ligand tests/test_data/ligand.pdbqt
```

## Project Structure

```
src/
  atom/           Atom types, van der Waals radii, H-bond properties
  molecule/       Molecule struct: atoms, bonds, torsions, geometry
  forcefield/     ForceField trait + Vina and AD4 implementations
  scoring/        High-level scoring combining forcefield terms
  optimization/   Optimizer trait + Monte Carlo search
  grid/           Spatial grid for efficient energy lookup
  io/             PDBQT parser and writer
  math/           Vector/matrix utilities
  utils/          General helpers
  lib.rs          Library root
  main.rs         CLI binary
tests/            Integration tests and test data
benches/          Criterion benchmarks
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
