# Docking parity over 5 seeds: rustdock-vina vs AutoDock Vina 1.2.7 (vina scoring)

CASF-2016 core set, subset `full` (285 complexes), docked with seeds 42, 1, 2, 3, 4; otherwise the protocol of the single-seed report. Rates are over all 285 complexes, as mean ± standard deviation over the seeds.

| | AutoDock Vina | rustdock-vina |
|---|---|---|
| Docked without error (%) | 98.2 ± 0.0 | 100.0 ± 0.0 |
| Top-1 success (< 2 A, %) | 53.3 ± 0.7 | 54.3 ± 0.2 |
| Any-mode success (< 2 A, %) | 78.2 ± 0.7 | 79.2 ± 1.7 |
| Median top-1 RMSD (A) | 1.65 ± 0.06 | 1.61 ± 0.04 |
| Top-1 PoseBusters-valid (%) | 90.7 ± 0.6 | 92.8 ± 0.4 |
| Median wall time (s) | 4.0 ± 0.0 | 3.0 ± 0.1 |
| Complexes whose top-1 outcome changes with the seed | 28 | 18 |

## Paired comparison (280 complexes both docked with every seed)

- Top-1 success, rustdock minus Vina: -0.1 percentage points (95% bootstrap CI -1.4 to +1.2; Wilcoxon signed-rank p = 0.37 on each complex's success rate over the seeds).
- Top-1 affinity, rustdock minus Vina: median -0.007 kcal/mol over all complexes and seeds; rustdock's is lower (better) in 21% and higher in 11% of docks (by more than 0.05 kcal/mol).

## Per seed

Top-1 successes, and McNemar's exact test on the complexes only one tool solved with that seed.

| Seed | Vina | rustdock | Only Vina | Only rustdock | McNemar p |
|---|---|---|---|---|---|
| 42 | 154/285 | 154/285 | 8 | 8 | 1.00 |
| 1 | 153/285 | 155/285 | 8 | 10 | 0.81 |
| 2 | 151/285 | 155/285 | 6 | 10 | 0.45 |
| 3 | 153/285 | 155/285 | 5 | 7 | 0.77 |
| 4 | 149/285 | 155/285 | 7 | 13 | 0.26 |
