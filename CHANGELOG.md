# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- Vina `xs` atom typing is now derived from the bond graph: carbons bonded to a
  heteroatom are `C_P` and not hydrophobic, and hydrogen-bond donors require an
  actual bonded polar hydrogen. Together these were worth about 5 kcal/mol of
  spurious attraction on the 1IEP benchmark.
- Ligand connectivity is now always inferred. Previously only the ~6 rotatable
  bonds from `BRANCH` records existed, so no bond-graph-dependent property
  could be computed.
- The ligand's intramolecular energy is now part of the search objective, so
  torsional moves can no longer fold the ligand through itself for free.
- Local optimization respects the search box; it could previously walk the
  ligand out of the box entirely and report the near-zero energy of a pose no
  longer touching the receptor.
- The L-BFGS optimizer was silently running as steepest descent (its curvature
  history could never populate) and its torsion degrees of freedom were inert.
  Both now work, and torsion angles no longer desynchronize from coordinates.
- Each search run now starts from a random pose in the box; box containment is
  tested on the ligand centroid rather than every atom, which no ligand larger
  than the box could ever satisfy.
- Binding modes are clustered by RMSD, so the reported modes are distinct
  rather than repeated copies of one pose.
- Unknown `--scoring` values, a missing `--ligand`, and malformed `--center`/
  `--size` are now errors instead of silent no-ops. Unrecognized config keys
  are reported.

### Changed
- Monte Carlo now minimizes each proposal before the Metropolis test, matching
  Vina's basin-hopping; temperature raised to Vina's 1.2.
- Logging defaults to `info`, and `dock` prints a Vina-style results table.
- `--scoring ad4` warns that it is experimental and unvalidated.

### Removed
- Placeholder `Molecule::calculate_internal_energy` (a bare inverse-distance
  sum) and the empty `optimization::genetic` module. `VinaScore` is now a
  working scorer rather than a stub.

## [0.1.0] - 2025-03-27

### Added
- Initial release
- Core molecular docking functionality
- Vina and (experimental) AutoDock4 scoring functions
- Parallelized search algorithm
- Command-line interface
- Test suite and documentation 