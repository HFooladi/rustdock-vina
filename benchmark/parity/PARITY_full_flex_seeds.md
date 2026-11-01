# Docking parity over 5 seeds: rustdock-vina vs AutoDock Vina 1.2.7 (vina scoring)

CASF-2016 core set, subset `full` (234 complexes), docked with seeds 42, 1, 2, 3, 4; otherwise the protocol of the single-seed report. Rates are over all 234 complexes, as mean ± standard deviation over the seeds.

Flexible side chains: up to 3 residues whose side chains have a heavy atom within 4 A of the crystal ligand, closest first, prepared with Meeko's `mk_prepare_receptor.py --flexres` (rigid part and flexible residues) from the same protein. The box covers the crystal ligand and the side chains plus 10 A. Complexes Meeko cannot prepare are left out for both programs (`flex_excluded.tsv` in the work directory). PoseBusters checks the top ligand pose against the protein with its crystal side chains.

| | AutoDock Vina | rustdock-vina |
|---|---|---|
| Docked without error (%) | 100.0 ± 0.0 | 100.0 ± 0.0 |
| Top-1 success (< 2 A, %) | 44.0 ± 1.2 | 44.4 ± 0.9 |
| Any-mode success (< 2 A, %) | 72.5 ± 1.3 | 72.9 ± 1.5 |
| Median top-1 RMSD (A) | 2.67 ± 0.12 | 2.63 ± 0.04 |
| Top-1 PoseBusters-valid (%) | 70.0 ± 0.8 | 69.5 ± 0.6 |
| Median wall time (s) | 15.6 ± 0.1 | 10.6 ± 0.1 |
| Complexes whose top-1 outcome changes with the seed | 32 | 25 |

## Paired comparison (234 complexes both docked with every seed)

- Top-1 success, rustdock minus Vina: +0.4 percentage points (95% bootstrap CI -1.1 to +2.1; Wilcoxon signed-rank p = 0.84 on each complex's success rate over the seeds).
- Top-1 affinity, rustdock minus Vina: median +0.000 kcal/mol over all complexes and seeds; rustdock's is lower (better) in 22% and higher in 22% of docks (by more than 0.05 kcal/mol).

## Per seed

Top-1 successes, and McNemar's exact test on the complexes only one tool solved with that seed.

| Seed | Vina | rustdock | Only Vina | Only rustdock | McNemar p |
|---|---|---|---|---|---|
| 42 | 101/234 | 104/234 | 3 | 6 | 0.51 |
| 1 | 106/234 | 106/234 | 10 | 10 | 1.00 |
| 2 | 100/234 | 101/234 | 5 | 6 | 1.00 |
| 3 | 106/234 | 103/234 | 10 | 7 | 0.63 |
| 4 | 102/234 | 106/234 | 4 | 8 | 0.39 |
