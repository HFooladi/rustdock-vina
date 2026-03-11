# Rustdock-Vina

A Rust implementation of [AutoDock Vina](https://vina.scripps.edu/) for molecular docking.

## Features

- Multiple scoring functions: Vina, AutoDock4 (AD4), and Vinardo
- Flexible side chain docking
- Parallelized Monte Carlo search via [rayon](https://github.com/rayon-rs/rayon)
- Configuration file support for docking parameters
- Batch ligand docking with directory output
- PDBQT file parsing and writing
- Memory-safe implementation in Rust

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
| `--flex <FILE>` | PDBQT file containing flexible receptor residues | |
| `-c, --config <FILE>` | Configuration file for docking parameters | |
| `--center <X,Y,Z>` | Center of the search box | |
| `--size <X,Y,Z>` | Size of the search box | |
| `--scoring <NAME>` | Scoring function: `vina`, `ad4`, or `vinardo` | `vina` |
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
| `--flex <FILE>` | PDBQT file containing flexible receptor residues | |
| `--scoring <NAME>` | Scoring function: `vina`, `ad4`, or `vinardo` | `vina` |

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
