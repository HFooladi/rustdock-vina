# Docking parity over 5 seeds: rustdock-vina vs AutoDock Vina 1.2.7 (vinardo scoring)

PoseBusters benchmark set (428 complexes), docked with seeds 42, 1, 2, 3, 4; otherwise the protocol of the single-seed report. Rates are over all 428 complexes, as mean ± standard deviation over the seeds.

| | AutoDock Vina | rustdock-vina |
|---|---|---|
| Docked without error (%) | 99.1 ± 0.0 | 100.0 ± 0.0 |
| Top-1 success (< 2 A, %) | 49.3 ± 0.6 | 49.5 ± 0.8 |
| Any-mode success (< 2 A, %) | 71.2 ± 0.9 | 70.6 ± 0.6 |
| Median top-1 RMSD (A) | 2.01 ± 0.02 | 2.04 ± 0.05 |
| Top-1 PoseBusters-valid (%) | 92.8 ± 0.7 | 93.1 ± 0.4 |
| Median wall time (s) | 5.7 ± 0.0 | 4.8 ± 0.0 |
| Complexes whose top-1 outcome changes with the seed | 55 | 59 |

## Paired comparison (424 complexes both docked with every seed)

- Top-1 success, rustdock minus Vina: -0.3 percentage points (95% bootstrap CI -1.5 to +0.8; Wilcoxon signed-rank p = 0.57 on each complex's success rate over the seeds).
- Top-1 affinity, rustdock minus Vina: median -0.001 kcal/mol over all complexes and seeds; rustdock's is lower (better) in 18% and higher in 15% of docks (by more than 0.05 kcal/mol).

## Per seed

Top-1 successes, and McNemar's exact test on the complexes only one tool solved with that seed.

| Seed | Vina | rustdock | Only Vina | Only rustdock | McNemar p |
|---|---|---|---|---|---|
| 42 | 210/428 | 215/428 | 13 | 18 | 0.47 |
| 1 | 214/428 | 213/428 | 16 | 15 | 1.00 |
| 2 | 211/428 | 206/428 | 20 | 15 | 0.50 |
| 3 | 208/428 | 211/428 | 15 | 18 | 0.73 |
| 4 | 213/428 | 214/428 | 18 | 19 | 1.00 |
