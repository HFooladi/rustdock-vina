# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- During Monte Carlo the ligand-receptor energy is capped at 10 per atom, as
  in Vina 1.2.7, not 1.5. 1.5 is the default of Vina's `monte_carlo`, but
  `Vina::global_search` sets `hunt_cap` to `(10, 10, 10)` before docking
  (`vina.cpp`). The parity tables are rerun with the fix.

### Added
- `scoring::score_pose`, which scores a pose with any force field and returns
  the breakdown `rustdock score` prints.
- A parity regression suite (`tests/parity.rs`) that runs in CI. It checks
  rustdock's crystal-pose scores (Vina and Vinardo) and its rescoring of
  Vina's docked poses against AutoDock Vina 1.2.7 to within 0.02 kcal/mol, on
  20 PoseBusters complexes. It also checks that docking 5 of them still finds
  the crystal pose. `benchmark/parity/make_golden.py` regenerates the fixtures.
- The parity harness docks with several seeds (`--seeds`) and several
  complexes at once on pinned cores (`--jobs`), and reports success rates
  over seeds with a paired comparison. The README's parity tables now come
  from five seeds; the Vina-scoring reports are rerun since the bond-inference
  fix.

## [0.2.0] - 2026-10-19

### Changed
- The command-line binary is now `rustdock` instead of `vina`, so it can be
  installed next to AutoDock Vina. Release archives are named
  `rustdock-<version>-<platform>`.
- `--flex` is now an error instead of being parsed and ignored; flexible
  receptor residues are not implemented.
- The docking search is now a port of AutoDock Vina 1.2.7's
  (`optimization::vina_search`). On the 29-complex CASF-2016 pilot, top-1
  success rose from 11/29 to 15/29 (Vina: 16/29), any-mode success from
  16/29 to 26/29 (Vina: 24/29), and the median top-1 RMSD fell from 3.94 to
  1.34 A (Vina: 1.31), at 0.72x Vina's total wall time. As in Vina:
  - `exhaustiveness` independent Monte Carlo chains, each running
    `70 * 3 * (50 + atoms + 10 * dof) / 2` steps, instead of
    `num_modes * exhaustiveness` short runs;
  - each step mutates one thing (position, orientation, or one torsion set
    to a random angle) and is minimized by Vina's dense BFGS and line search
    for `(25 + atoms) / 3` steps before the Metropolis test at 1.2;
  - grids span exactly the search box, and atoms outside it are clamped and
    penalized; clashes are softened during the search (`hunt_cap`);
  - each chain keeps its distinct low minima; chains merge at 2 A, modes
    are deduplicated at `min_rmsd`, then refined and rescored on exact
    energies;
  - a mode's energy is its intermolecular energy plus its internal energy
    relative to the best mode's.
- Intramolecular pair energies are tabulated per pair of atom classes
  (as Vina's `precalculate`); they were 90% of each evaluation.
- `MonteCarloParams` now holds Vina's settings (`exhaustiveness`,
  `temperature`, `mutation_amplitude`, `hunt_cap`, `min_rmsd`,
  `local_steps`); the old per-move probabilities and step limits are gone.
  In `Optimizer::generate_poses`, `num_poses` is the number of modes and
  `max_steps` the Monte Carlo steps per chain (`0` for Vina's default).
- Monte Carlo now minimizes each proposal before the Metropolis test, matching
  Vina's basin-hopping; temperature raised to Vina's 1.2.
- Logging defaults to `info`, and `dock` prints a Vina-style results table.
- `--scoring ad4` warns that it is experimental and unvalidated.
- The local optimizer rotates about the entry pose's centroid rather than the
  current one, so rotation and torsion derivatives decouple.
- Affinity grids are built ~5x faster (1IEP: 169 ms to 31 ms; a 1-run dock
  0.5 s to 0.2 s). Pair energies come from radial tables per receptor/ligand
  class, as in Vina's precalculation, and all ligand classes are filled in one
  sweep that only visits receptor atoms within the cutoff. Grid values match
  exact scoring to ~2e-4 relative.

### Added
- Batch screening: `dock --batch <files or dirs> --dir <out>` docks many
  ligands into one box in parallel, sharing one set of receptor grids, and
  writes `scores.tsv` with every best affinity. A ligand that fails is
  listed with its error and skipped. `--cpu` sets the thread count. Each
  ligand's poses are identical to docking it alone with the same seed. 100
  CASF-2016 ligands into 1e66's pocket on 96 threads took 141 s, against
  1024 s for Vina 1.2.7's `--batch` and 210 s for 12 parallel Vina
  processes of 8 threads; best affinities agree with Vina's at Pearson
  0.98. Library API: `MonteCarlo::target` and `MonteCarlo::dock_into`, with
  grids shared through `grid::ReceptorGrids`.
- Vinardo scoring (`--scoring vinardo`, `VinardoForceField`), following
  AutoDock Vina 1.2.7: Vina's typing and torsional term with Vinardo's radii,
  one Gaussian, and its hydrophobic and H-bond ramps. On the CASF-2016 pilot
  it rescores crystal poses within 0.02 kcal/mol of Vina's Vinardo and docks
  19/29 top-1 (Vina: 18/29). The parity harness takes `--scoring vinardo`.
  On the full sets, top-1 success is CASF-2016 153/285 (Vina: 144/285) and
  PoseBusters 222/428 (Vina: 210/428), crystal poses rescore within 0.03
  kcal/mol of Vina's, and Vina rescores rustdock's output within 0.011
  (`PARITY_full_vinardo.md`, `PARITY_posebusters_vinardo.md`).
- `benchmark/parity/`: a CASF-2016 docking parity harness against AutoDock
  Vina 1.2.7. Both programs dock identical Meeko/Open Babel-prepared inputs;
  the report gives top-1 success, symmetry-aware RMSD, PoseBusters validity,
  wall time and crystal-pose rescoring per complex. Tools install with uv.
  It also runs the PoseBusters benchmark set (`--subset posebusters`). Top-1
  success: CASF-2016 154/285 (Vina: 154/285), PoseBusters 238/428 (Vina:
  236/428), at 0.7x and 0.8x Vina's median wall time
  (`PARITY_full.md`, `PARITY_posebusters.md`).
- `dock --seed` makes a docking run reproducible.
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
- Bonds are inferred as in Vina's `model::assign_bonds`: closer than 1.1
  times the sum of the covalent radii (not a flat 2.1 A, or 1.3 A to
  hydrogen), and not when a third atom is closer to both. In 2wnc the amine
  hydrogen sits between the N and its methyl carbon, so Vina types that
  carbon hydrophobic. Crystal-pose rescoring of the CASF-2016 and
  PoseBusters sets went from max 0.38 to 0.023 kcal/mol off Vina's.
- Macrocycles opened by Meeko now dock and score as in Vina 1.2.7. The
  closure types `CG0`-`CG3` and `G0`-`G3` are read (the atom type is read to
  the end of the line, not two columns); bonds are inferred only within a
  rigid fragment, so a dummy no longer bonds to the carbon it sits on; and
  Vina's closure rules apply: dummies are not scored against the receptor or
  other atoms, pairs within three bonds of the closure are excluded, and
  during the search each dummy is pulled onto its carbon by a `50 * r`
  attraction. Crystal-pose scores of the 12 CASF-2016 macrocycles went from
  up to 1.23 kcal/mol off Vina's to within 0.013.
- During Monte Carlo the ligand-receptor energy is capped at 1.5 per atom
  (`MonteCarloParams::hunt_cap_receptor`, Vina's `hunt_cap[1]`), not 10.
- PDBQT output now round-trips. Docking results are written in Vina's
  layout (`MODEL n`, `REMARK VINA RESULT: affinity rmsd_lb rmsd_ub`, then
  `INTER + INTRA`, `INTER`, `INTRA`, `UNBOUND`) followed by the ligand's own
  records with only the coordinates changed, so Meeko's `REMARK SMILES`
  lines, the `ROOT`/`BRANCH` tree, `TORSDOF`, atom names and charges survive.
  Previously all of these were dropped and the atom columns were shifted. A
  Meeko ligand written back unchanged is byte-identical; on the CASF pilot,
  Meeko reads every rustdock output and Vina rescores rustdock's top pose
  within 0.007 kcal/mol of the affinity rustdock reports. New
  `Molecule::pdbqt_records` holds the records.
- Vina scoring now matches AutoDock Vina 1.2.7 exactly. Rescoring the 29
  CASF-2016 pilot crystal poses, the median gap fell from 1.73 to 0.003
  kcal/mol (max 6.57 to 0.012), and pilot top-1 docking success rose from
  8/29 to 11/29 (Vina: 16/29). Typing now follows Vina's
  `model::assign_types`:
  - any N or O bonded to a polar hydrogen is a donor, not only `NA`/`OA`,
    which had dropped every backbone N-H and amide/amine donor;
  - metals are donors;
  - `SA` is not an acceptor;
  - halogens and metals count as heteroatoms, so carbons bonded to them are
    polar, not hydrophobic;
  - sodium (`Na`) parsed as an acceptor nitrogen (`NA`); type names are now
    case-sensitive like Vina's, which adds Vina's other metals (K, Cu, Ni, Co,
    Cd, Hg, U).
- Affinity divides by Vina's weighted torsion count
  (`Molecule::vina_torsion_count`: torsions that only spin a terminal group
  count 0.5) instead of the number of torsions. `calculate_affinity` and
  `ForceField::apply_torsional_normalization` take it as `f64`. `dock` and
  `score` use the force field's normalization rather than Vina's weight for
  every force field.
- PDBQT files with negative residue numbers (legal in PDB, and written by
  Open Babel) failed to parse. `Atom::residue_num` is now `i32`.
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