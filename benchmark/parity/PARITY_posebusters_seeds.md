# Docking parity over 5 seeds: rustdock-vina vs AutoDock Vina 1.2.7 (vina scoring)

PoseBusters benchmark set (428 complexes), docked with seeds 42, 1, 2, 3, 4; otherwise the protocol of the single-seed report. Rates are over all 428 complexes, as mean ± standard deviation over the seeds.

| | AutoDock Vina | rustdock-vina |
|---|---|---|
| Docked without error (%) | 99.1 ± 0.0 | 100.0 ± 0.0 |
| Top-1 success (< 2 A, %) | 54.7 ± 0.8 | 55.2 ± 0.3 |
| Any-mode success (< 2 A, %) | 77.6 ± 1.4 | 79.0 ± 0.8 |
| Median top-1 RMSD (A) | 1.77 ± 0.03 | 1.77 ± 0.01 |
| Top-1 PoseBusters-valid (%) | 93.4 ± 0.5 | 93.7 ± 0.6 |
| Median wall time (s) | 5.3 ± 0.3 | 4.5 ± 0.1 |
| Complexes whose top-1 outcome changes with the seed | 58 | 45 |

## Paired comparison (424 complexes both docked with every seed)

- Top-1 success, rustdock minus Vina: -0.2 percentage points (95% bootstrap CI -1.6 to +1.1; Wilcoxon signed-rank p = 0.71 on each complex's success rate over the seeds).
- Top-1 affinity, rustdock minus Vina: median -0.006 kcal/mol over all complexes and seeds; rustdock's is lower (better) in 23% and higher in 12% of docks (by more than 0.05 kcal/mol).

## Per seed

Top-1 successes, and McNemar's exact test on the complexes only one tool solved with that seed.

| Seed | Vina | rustdock | Only Vina | Only rustdock | McNemar p |
|---|---|---|---|---|---|
| 42 | 236/428 | 236/428 | 17 | 17 | 1.00 |
| 1 | 237/428 | 235/428 | 17 | 15 | 0.86 |
| 2 | 230/428 | 237/428 | 12 | 19 | 0.28 |
| 3 | 237/428 | 238/428 | 14 | 15 | 1.00 |
| 4 | 231/428 | 235/428 | 17 | 21 | 0.63 |
