# Virtual-screening parity (DUDE-Z)

`screen.py` asks whether rustdock ranks actives above decoys as well as
AutoDock Vina 1.2.7 does. The pose benchmarks in `../parity` check that the
top pose is right; this benchmark checks the scores across thousands of
ligands, which is how a virtual screen uses them.

## Data

[DUDE-Z](https://dudez.docking.org) (Stein et al., *J. Chem. Inf. Model.*
2021, [10.1021/acs.jcim.0c00598](https://doi.org/10.1021/acs.jcim.0c00598)),
the successor of DUD-E. It is downloaded on first use, for each target:

- `DOCKING_GRIDS_AND_POSES/<T>/rec.crg.pdb`, the prepared receptor with its
  polar hydrogens, and `xtal-lig.pdb`, its crystal ligand;
- `DUDE-Z-benchmark-grids/<T>/ligands.smi` (actives from ChEMBL, often in
  several protomers) and `decoys.smi` (50 property-matched ZINC decoys per
  active).

The default is ten targets of different families: AA2AR, ABL1, ACES, AMPC,
ANDR, EGFR, HIVPR, HS90A, PARP1 and THRB.

## Protocol

- Receptor: `rec.crg.pdb` to PDBQT with Open Babel (`-xr`).
- Box: the crystal ligand's heavy-atom extent plus 10 Å per axis, centred on
  it, as in `../parity`.
- Ligands: RDKit ETKDG conformer (seed 42) and MMFF, then Meeko, from the
  SMILES as given (DUDE-Z's protonation). Ligands that fail are left out.
- Docking: both programs get the same files and settings: Vina scoring,
  exhaustiveness 8, nine modes, one thread per dock with each worker pinned to
  its own physical core, and docking seeds 42 and 1. Jobs of both programs
  are shuffled together so they share the load.
- A compound scores the best top-pose affinity over its protomers. Only
  compounds that every run docked are ranked.

## Metrics

Per target: ROC AUC, BEDROC (α = 20; Truchon and Bayly 2007) and enrichment
factors at 1% and 5%. Ties are ranked decoy first, so they never help.

The report compares rustdock with Vina at seed 42, and sets the difference
against each program's own seed-to-seed difference (seed 1 vs seed 42). The
docking search is stochastic, so two Vina runs disagree too; a rustdock − Vina
difference within that range is noise. Confidence intervals come from
resampling actives and decoys within each target (2,000 draws).

## Run

From the repository root, with the venv and Vina binary of `../parity`:

```bash
P=benchmark/parity/.venv/bin/python
$P -I benchmark/screening/screen.py run --work /path/to/screening --jobs 64
# pilot: two targets, eight actives and eight decoys each
$P -I benchmark/screening/screen.py run --work /tmp/pilot --targets AA2AR AMPC --limit 8
# report again from finished docks
$P -I benchmark/screening/screen.py report --work /path/to/screening --report SCREENING.md
```

Finished docks are cached as `<work>/<target>/docked/*.json`, so a stopped run
resumes. The full set (41,252 ligands, four runs, 165,008 docks) took 18 hours on
64 cores of an AMD EPYC 7713; `--decoys N` samples N decoys per target instead.

## Results (2026-11-04)

Ten targets, 665 actives and 40,443 decoys, rustdock at 57a500e (v0.3.0).
Full report: [`SCREENING.md`](SCREENING.md).

| Mean over 10 targets | Vina | rustdock | rustdock − Vina (95% CI) | Vina − Vina, two seeds (95% CI) |
|---|---|---|---|---|
| ROC AUC | 0.678 | 0.678 | +0.000 (−0.003 to +0.003) | +0.002 (−0.001 to +0.005) |
| BEDROC (α = 20) | 0.192 | 0.189 | −0.003 (−0.008 to +0.001) | −0.001 (−0.006 to +0.005) |
| EF1% | 7.1 | 6.7 | −0.37 (−1.21 to +0.69) | −0.22 (−0.95 to +0.95) |
| EF5% | 3.76 | 3.58 | −0.19 (−0.38 to +0.09) | −0.09 (−0.30 to +0.16) |
| Mean seconds per dock (one thread) | 26.7 | 23.8 | | |

- **Ranking.** Every difference is inside its CI around zero and of the size of
  Vina's own seed-to-seed difference. Averaged over both seeds, rustdock has
  the higher AUC on 5 targets and Vina on 5; for EF5%, Vina is ahead on 7 of
  10 (Wilcoxon p = 0.07), by 0.18 on 3.7, which the CIs do not separate from
  noise.
- **Scores.** Per compound, rustdock's best affinity correlates with Vina's
  at Spearman 0.994, the same as Vina with itself across seeds (0.994), on
  every target. 98.3% of compounds are within 0.5 kcal/mol of Vina's score
  (Vina vs Vina: 98.4%).
- **Failures.** No dock failed in either program; 7 of the 41,259 SMILES
  could not be embedded by RDKit and are left out of both.
- **Speed.** rustdock was 11% faster per dock on one thread, and more so for
  flexible ligands (113 s vs 140 s above 10 torsions).

## rustdock finds slightly lower minima

rustdock's top affinity is more than 0.05 kcal/mol below Vina's for 19.9% of
compounds, and above it for 11.1%; two Vina seeds split 13.7/13.7. The mean
shift is −0.016 kcal/mol, alike for actives (−0.012) and decoys (−0.016), so
it moves the whole ranking and does not change it.

- **It is real under Vina's scoring.** For 500 top poses, Vina's
  `--score_only` gives rustdock's reported affinity to within 0.004 kcal/mol
  on average (0.023 at most), and 94% of the poses where rustdock reported a
  lower energy are lower under Vina's scoring too.
- **It is mostly local.** Both programs put the top pose within 1 Å of each
  other as often as two Vina seeds do (75%). Within the same pose, rustdock's
  energy is lower in 47% of compounds and Vina's in 25%, by a median 0.008
  kcal/mol; two Vina seeds split 33/33. It appears at every flexibility,
  from 0–2 torsions (14.7% vs 5.0%) to more than 10.
- **The cause is Vina's lookup tables.** Vina 1.2.7 evaluates pair energies
  from tables sampled 32 times per Å², interpolating the energy and its
  derivative separately, also in the final exact-energy refinement
  (`non_cache`). rustdock computes them exactly there. Re-docking 2,068
  compounds with Vina's tables emulated:

| rustdock variant | same pose: rustdock lower | same pose: Vina lower | > 0.05 kcal/mol: rustdock lower | > 0.05 kcal/mol: Vina lower |
|---|---|---|---|---|
| Stock | 46.9% | 23.2% | 20.9% | 10.5% |
| Vina's ligand-internal tables | 45.4% | 25.6% | 20.6% | 11.4% |
| Vina's grid node values | 46.5% | 25.6% | 20.4% | 10.8% |
| Vina's tables in the final refinement | 40.7% | 27.3% | 16.6% | 12.4% |
| All three | 40.8% | 29.0% | 14.9% | 13.7% |
| *Vina seed 1 vs seed 42 (noise)* | *31.6%* | *34.1%* | *13.3%* | *13.3%* |

  With all three, the split of larger differences is that of Vina against
  itself, and the same-pose median gap halves (−0.004 kcal/mol). rustdock keeps
  its exact energies: they are the better approximation of the same scoring
  function, and they change neither the ranking nor the pose.
