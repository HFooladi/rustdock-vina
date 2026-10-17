# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- The docking search now scores ligand-receptor energy from precomputed
  per-atom-class affinity grids (`GridSet`), as Vina does, instead of summing
  over receptor atoms for every pose. The final polish of each run and all
  reported energies still use exact pairwise scoring. On the 1IEP redocking
  benchmark a run is ~6x cheaper, and at equal wall-clock the extra runs raise
  top-pose recovery (< 2 A) from 2/16 to 7/16 seeds. New `dock` flags:
  `--grid-spacing` (default 0.375 A) and `--no-grid`. Force fields opt in via
  `ForceField::grid_cutoff`; AD4 (charge-dependent) stays pairwise.
- The local optimizer computes its gradient analytically (per-atom forces
  folded into force, torque and per-torsion projections, as Vina does) instead
  of by finite differences, which cost `7 + n_torsions` energy evaluations per
  gradient. New `ForceField::atom_pair_derivative` (closed form for Vina,
  central difference by default). On 1IEP a 64-run grid dock drops from 1.4 s
  to 0.6 s and an exact dock from 8.0 s to 2.0 s; at 128 runs (0.7 s) the top
  pose is within 2 A of the crystal pose for 15/16 seeds.

### Fixed
- Per-run seeds in `generate_poses` were `master + i`, so docks with nearby
  master seeds repeated most of their runs (masters 0-15 over 8 runs shared
  just 23 distinct seeds). Run seeds are now mixed with splitmix64.
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
- The local optimizer rotates about the entry pose's centroid rather than the
  current one, so rotation and torsion derivatives decouple.

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