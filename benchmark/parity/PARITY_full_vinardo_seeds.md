# Docking parity over 5 seeds: rustdock-vina vs AutoDock Vina 1.2.7 (vinardo scoring)

CASF-2016 core set, subset `full` (285 complexes), docked with seeds 42, 1, 2, 3, 4; otherwise the protocol of the single-seed report. Rates are over all 285 complexes, as mean ± standard deviation over the seeds.

| | AutoDock Vina | rustdock-vina |
|---|---|---|
| Docked without error (%) | 98.2 ± 0.0 | 100.0 ± 0.0 |
| Top-1 success (< 2 A, %) | 51.5 ± 1.1 | 52.9 ± 1.2 |
| Any-mode success (< 2 A, %) | 72.2 ± 2.8 | 73.5 ± 0.9 |
| Median top-1 RMSD (A) | 1.84 ± 0.07 | 1.82 ± 0.07 |
| Top-1 PoseBusters-valid (%) | 90.8 ± 0.8 | 92.0 ± 0.6 |
| Median wall time (s) | 4.3 ± 0.1 | 3.6 ± 0.0 |
| Complexes whose top-1 outcome changes with the seed | 39 | 34 |

## Paired comparison (280 complexes both docked with every seed)

- Top-1 success, rustdock minus Vina: +0.8 percentage points (95% bootstrap CI -0.5 to +2.1; Wilcoxon signed-rank p = 0.28 on each complex's success rate over the seeds).
- Top-1 affinity, rustdock minus Vina: median -0.002 kcal/mol over all complexes and seeds; rustdock's is lower (better) in 15% and higher in 14% of docks (by more than 0.05 kcal/mol).

## Per seed

Top-1 successes, and McNemar's exact test on the complexes only one tool solved with that seed.

| Seed | Vina | rustdock | Only Vina | Only rustdock | McNemar p |
|---|---|---|---|---|---|
| 42 | 144/285 | 155/285 | 4 | 15 | 0.02 |
| 1 | 146/285 | 146/285 | 9 | 9 | 1.00 |
| 2 | 152/285 | 149/285 | 10 | 7 | 0.63 |
| 3 | 146/285 | 153/285 | 6 | 13 | 0.17 |
| 4 | 146/285 | 151/285 | 11 | 16 | 0.44 |
