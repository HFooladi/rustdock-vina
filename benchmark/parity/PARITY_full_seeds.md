# Docking parity over 5 seeds: rustdock-vina vs AutoDock Vina 1.2.7 (vina scoring)

CASF-2016 core set, subset `full` (285 complexes), docked with seeds 42, 1, 2, 3, 4; otherwise the protocol of the single-seed report. Rates are over all 285 complexes, as mean ± standard deviation over the seeds.

| | AutoDock Vina | rustdock-vina |
|---|---|---|
| Docked without error (%) | 98.2 ± 0.0 | 100.0 ± 0.0 |
| Top-1 success (< 2 A, %) | 53.3 ± 0.7 | 53.3 ± 0.8 |
| Any-mode success (< 2 A, %) | 78.2 ± 0.7 | 77.6 ± 1.6 |
| Median top-1 RMSD (A) | 1.65 ± 0.06 | 1.74 ± 0.05 |
| Top-1 PoseBusters-valid (%) | 90.7 ± 0.6 | 92.6 ± 0.8 |
| Median wall time (s) | 4.0 ± 0.0 | 3.1 ± 0.0 |
| Complexes whose top-1 outcome changes with the seed | 28 | 19 |

## Paired comparison (280 complexes both docked with every seed)

- Top-1 success, rustdock minus Vina: -1.1 percentage points (95% bootstrap CI -2.6 to +0.3; Wilcoxon signed-rank p = 0.02 on each complex's success rate over the seeds).
- Top-1 affinity, rustdock minus Vina: median -0.007 kcal/mol over all complexes and seeds; rustdock's is lower (better) in 18% and higher in 11% of docks (by more than 0.05 kcal/mol).

## Per seed

Top-1 successes, and McNemar's exact test on the complexes only one tool solved with that seed.

| Seed | Vina | rustdock | Only Vina | Only rustdock | McNemar p |
|---|---|---|---|---|---|
| 42 | 154/285 | 150/285 | 9 | 5 | 0.42 |
| 1 | 153/285 | 153/285 | 6 | 6 | 1.00 |
| 2 | 151/285 | 155/285 | 6 | 10 | 0.45 |
| 3 | 153/285 | 152/285 | 9 | 8 | 1.00 |
| 4 | 149/285 | 149/285 | 9 | 9 | 1.00 |
