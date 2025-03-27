# Rustdock-Vina

A Rust implementation of AutoDock Vina for molecular docking.

## Overview

Rustdock-Vina is an open-source reimplementation of the popular AutoDock Vina molecular docking software in Rust. It aims to provide improved performance, memory safety, and a more maintainable codebase while preserving the original algorithm's accuracy.

Key features:
- Modern, memory-safe implementation in Rust
- Multiple scoring functions (Vina, AutoDock4)
- Flexible side chain docking
- Parallelized search algorithm for better performance
- User-friendly command-line interface

## Installation

### From Cargo (Rust package manager)


## Installation

### Prerequisites

- Rust 1.70 or later
- Cargo (Rust's package manager)

### Building from source

```bash
git clone https://github.com/yourusername/rustdock-vina
cd rustdock-vina
cargo build --release
```

The compiled binary will be available at `target/release/rustdock-vina`.

## Usage

Basic usage example:

```bash
rustdock-vina --receptor protein.pdbqt --ligand ligand.pdbqt --center 0,0,0 --size 20,20,20
```

For more detailed usage instructions, please refer to the [documentation](docs/).

## Documentation

- [User Guide](docs/user-guide.md)
- [API Documentation](docs/api.md)
- [Contributing Guidelines](CONTRIBUTING.md)

## Testing

Run the test suite:

```bash
cargo test
```

Run benchmarks:

```bash
cargo bench
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please read our [Contributing Guidelines](CONTRIBUTING.md) for details on how to submit pull requests, report issues, and contribute to the project.

## Acknowledgments

- AutoDock Vina developers for the original implementation and research
- The Rust community for excellent tools and libraries 