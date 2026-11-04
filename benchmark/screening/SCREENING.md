# Virtual-screening parity: rustdock vs AutoDock Vina 1.2.7

DUDE-Z, 10 targets, 665 actives and 40443 decoys (all per target), docked by both programs with the same Meeko PDBQTs, box (crystal ligand + 10 A) and settings: Vina scoring, exhaustiveness 8, one thread per dock on its own core, seeds 42, 1. A compound scores the best top-pose affinity over its protomers; compounds that failed in any run are left out of all.

## Per target (seed 42)

| Target | Actives | Decoys | AUC Vina | AUC rustdock | BEDROC Vina | BEDROC rustdock | EF1% Vina | EF1% rustdock | Spearman ρ |
|---|---|---|---|---|---|---|---|---|---|
| AA2AR | 87 | 4500 | 0.684 | 0.681 | 0.158 | 0.149 | 3.4 | 2.3 | 0.996 |
| ABL1 | 49 | 3150 | 0.727 | 0.725 | 0.200 | 0.199 | 8.2 | 8.2 | 0.996 |
| ACES | 84 | 4949 | 0.616 | 0.623 | 0.112 | 0.109 | 1.2 | 0.0 | 0.991 |
| AMPC | 48 | 3148 | 0.593 | 0.592 | 0.080 | 0.071 | 0.0 | 2.1 | 0.991 |
| ANDR | 94 | 4950 | 0.616 | 0.616 | 0.301 | 0.303 | 17.2 | 17.2 | 0.999 |
| EGFR | 68 | 4248 | 0.698 | 0.702 | 0.129 | 0.128 | 4.4 | 4.4 | 0.994 |
| HIVPR | 46 | 3200 | 0.581 | 0.581 | 0.088 | 0.072 | 2.2 | 0.0 | 0.993 |
| HS90A | 50 | 3500 | 0.563 | 0.557 | 0.029 | 0.030 | 0.0 | 0.0 | 0.981 |
| PARP1 | 95 | 6149 | 0.873 | 0.879 | 0.478 | 0.487 | 25.4 | 26.5 | 0.994 |
| THRB | 44 | 2649 | 0.831 | 0.828 | 0.348 | 0.345 | 9.1 | 6.8 | 0.976 |
| **Mean** | | | 0.678 | 0.678 | 0.192 | 0.189 | 7.1 | 6.7 | |

## rustdock − Vina, against Vina's seed-to-seed noise

Mean over targets of the paired difference, 95% CI from resampling compounds within each target. A rustdock − Vina difference inside the range of Vina − Vina (same program, another seed) is noise.

| Metric | rustdock − Vina | Vina − Vina (seeds) | rustdock − rustdock (seeds) |
|---|---|---|---|
| AUC | +0.000 (-0.003 to +0.003) | +0.002 (-0.001 to +0.005) | -0.001 (-0.005 to +0.002) |
| BEDROC | -0.003 (-0.008 to +0.001) | -0.001 (-0.006 to +0.005) | -0.002 (-0.007 to +0.004) |
| EF1% | -0.368 (-1.212 to +0.689) | -0.222 (-0.949 to +0.946) | +0.161 (-0.723 to +0.977) |
| EF5% | -0.186 (-0.379 to +0.092) | -0.085 (-0.297 to +0.157) | -0.083 (-0.278 to +0.121) |

## Per-compound scores

| Pair | Spearman ρ | Mean diff (kcal/mol) | Mean abs diff | Within 0.5 kcal/mol |
|---|---|---|---|---|
| rustdock.s42 − vina.s42 | 0.994 | -0.016 | 0.070 | 98.3% |
| vina.s1 − vina.s42 | 0.994 | -0.002 | 0.065 | 98.4% |
| rustdock.s1 − rustdock.s42 | 0.995 | -0.000 | 0.051 | 98.6% |

## Time and failures

| Run | Mean seconds per dock (1 thread) | Failed docks |
|---|---|---|
| vina.s42 | 26.7 | 0 |
| vina.s1 | 26.7 | 0 |
| rustdock.s42 | 23.8 | 0 |
| rustdock.s1 | 23.8 | 0 |
