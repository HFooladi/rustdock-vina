# Docking parity over 5 seeds: rustdock-vina vs AutoDock Vina 1.2.7 (vinardo scoring)

CASF-2016 core set, subset `full` (285 complexes), docked with seeds 42, 1, 2, 3, 4; otherwise the protocol of the single-seed report. Rates are over all 285 complexes, as mean ± standard deviation over the seeds.

| | AutoDock Vina | rustdock-vina |
|---|---|---|
| Docked without error (%) | 98.2 ± 0.0 | 100.0 ± 0.0 |
| Top-1 success (< 2 A, %) | 51.5 ± 1.1 | 54.0 ± 0.8 |
| Any-mode success (< 2 A, %) | 72.2 ± 2.8 | 75.3 ± 1.3 |
| Median top-1 RMSD (A) | 1.84 ± 0.07 | 1.76 ± 0.05 |
| Top-1 PoseBusters-valid (%) | 90.8 ± 0.8 | 92.4 ± 0.8 |
| Median wall time (s) | 4.3 ± 0.1 | 3.3 ± 0.0 |
| Complexes whose top-1 outcome changes with the seed | 39 | 21 |

## Paired comparison (280 complexes both docked with every seed)

- Top-1 success, rustdock minus Vina: +2.2 percentage points (95% bootstrap CI +0.6 to +3.9; Wilcoxon signed-rank p = 0.02 on each complex's success rate over the seeds).
- Top-1 affinity, rustdock minus Vina: median -0.003 kcal/mol over all complexes and seeds; rustdock's is lower (better) in 18% and higher in 9% of docks (by more than 0.05 kcal/mol).

## Per seed

Top-1 successes, and McNemar's exact test on the complexes only one tool solved with that seed.

| Seed | Vina | rustdock | Only Vina | Only rustdock | McNemar p |
|---|---|---|---|---|---|
| 42 | 144/285 | 153/285 | 7 | 16 | 0.09 |
| 1 | 146/285 | 157/285 | 5 | 16 | 0.03 |
| 2 | 152/285 | 151/285 | 7 | 6 | 1.00 |
| 3 | 146/285 | 155/285 | 6 | 15 | 0.08 |
| 4 | 146/285 | 154/285 | 6 | 14 | 0.12 |
