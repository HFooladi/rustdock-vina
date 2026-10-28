# Docking parity over 5 seeds: rustdock-vina vs AutoDock Vina 1.2.7 (vinardo scoring)

PoseBusters benchmark set (428 complexes), docked with seeds 42, 1, 2, 3, 4; otherwise the protocol of the single-seed report. Rates are over all 428 complexes, as mean ± standard deviation over the seeds.

| | AutoDock Vina | rustdock-vina |
|---|---|---|
| Docked without error (%) | 99.1 ± 0.0 | 100.0 ± 0.0 |
| Top-1 success (< 2 A, %) | 49.3 ± 0.6 | 51.6 ± 0.4 |
| Any-mode success (< 2 A, %) | 71.2 ± 0.9 | 74.3 ± 0.8 |
| Median top-1 RMSD (A) | 2.01 ± 0.02 | 1.91 ± 0.02 |
| Top-1 PoseBusters-valid (%) | 92.8 ± 0.7 | 93.6 ± 0.5 |
| Median wall time (s) | 5.7 ± 0.0 | 4.6 ± 0.0 |
| Complexes whose top-1 outcome changes with the seed | 55 | 50 |

## Paired comparison (424 complexes both docked with every seed)

- Top-1 success, rustdock minus Vina: +1.8 percentage points (95% bootstrap CI +0.5 to +3.1; Wilcoxon signed-rank p = 0.00 on each complex's success rate over the seeds).
- Top-1 affinity, rustdock minus Vina: median -0.002 kcal/mol over all complexes and seeds; rustdock's is lower (better) in 20% and higher in 12% of docks (by more than 0.05 kcal/mol).

## Per seed

Top-1 successes, and McNemar's exact test on the complexes only one tool solved with that seed.

| Seed | Vina | rustdock | Only Vina | Only rustdock | McNemar p |
|---|---|---|---|---|---|
| 42 | 210/428 | 222/428 | 8 | 20 | 0.04 |
| 1 | 214/428 | 222/428 | 10 | 18 | 0.18 |
| 2 | 211/428 | 220/428 | 11 | 20 | 0.15 |
| 3 | 208/428 | 222/428 | 7 | 21 | 0.01 |
| 4 | 213/428 | 218/428 | 14 | 19 | 0.49 |
