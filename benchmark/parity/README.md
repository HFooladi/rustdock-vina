# CASF-2016 docking parity

`parity.py` docks the same complexes with rustdock-vina and AutoDock Vina 1.2.7
from identical inputs, and writes a comparison table (`PARITY_<subset>.md`).

## Setup

Python tools go in a uv venv; Vina is the official release binary.

```bash
cd benchmark/parity
uv venv --python 3.12 .venv
uv pip install --python .venv -r requirements.txt
gh release download v1.2.7 -R ccsb-scripps/AutoDock-Vina -p vina_1.2.7_linux_x86_64 -D bin
mv bin/vina_1.2.7_linux_x86_64 bin/vina && chmod +x bin/vina
```

Open Babel 3.x (`obabel`) must be on `PATH` for receptor preparation.

Structures come from a local PDBbind copy (any version containing the CASF-2016
core set), laid out as `<root>/<group>/<code>/<code>_{protein.pdb,ligand.sdf}`.

The PoseBusters benchmark set is `posebusters_benchmark_set/` in
`posebusters_paper_data.zip` from [Zenodo record 8278563](https://zenodo.org/records/8278563)
(54 MB). Its proteins have no hydrogens, so the harness adds them with Open
Babel at pH 7.4 for both programs.

## Run

```bash
benchmark/parity/.venv/bin/python benchmark/parity/parity.py run \
    --pdbbind /path/to/PDBbind --subset pilot --cpu 16
```

Subsets of `casf2016_coreset.tsv`:

| Subset | Complexes | Selection |
|---|---|---|
| `pilot` | 29 | Best-resolved complex of each odd-numbered target |
| `targets` | 57 | Best-resolved complex of every target |
| `full` | 285 | The whole core set |
| `posebusters` | 428 | Every complex in the PoseBusters benchmark set (pass its directory as `--data`) |

Steps cache under `work/<code>/`, so a rerun resumes where it stopped.
`parity.py report` rebuilds the table from `work/results.csv`.

`--scoring vinardo` compares both programs' Vinardo instead (work directory
`work-vinardo/`, report `PARITY_<subset>_vinardo.md`).

### Several seeds

One seed cannot tell a real difference of a few complexes from chance. To
dock every complex with several seeds, use `--seeds`. To dock several
complexes at once, use `--jobs`; each job is pinned to its own `--cpu` cores.

```bash
benchmark/parity/.venv/bin/python benchmark/parity/parity.py run \
    --pdbbind /path/to/PDBbind --subset full --seeds 42 1 2 3 4 --jobs 4 --cpu 16
```

Seed 42 keeps the file names above. Other seeds add `_s<seed>` to the docking
outputs and to `results_s<seed>.csv`. The run also writes
`PARITY_<subset>_seeds[_vinardo].md`, which `parity.py seeds --seeds ...`
rebuilds from the CSVs. That report gives:

- each rate as mean ± standard deviation over the seeds;
- a paired comparison over complexes: each complex's top-1 success rate across
  the seeds, with a bootstrap 95% CI of the difference and a Wilcoxon
  signed-rank test;
- McNemar's exact test for each seed.

With `--jobs`, choose the slots so that no two jobs share a physical core;
otherwise wall times are inflated. On a 64-core machine with 2-way SMT, use
`--jobs 4 --cpu 16`.

## Protocol

| Step | How |
|---|---|
| Receptor | `protein.pdb` without waters → Open Babel `-xr` PDBQT (hydrogens added at pH 7.4 when the file has none) |
| Ligand | Crystal SDF → RDKit ETKDG conformer (seed 42) → Meeko PDBQT, so the crystal pose never reaches the search |
| Box | Crystal ligand centroid; size = heavy-atom extent + 10 Å per axis |
| Search | Both tools: exhaustiveness 8, 9 modes, seed 42, same CPU count, their own defaults otherwise |
| Pose RMSD | Symmetry-aware heavy-atom RMSD to the crystal ligand, in place (RDKit `CalcRMS`) |
| Validity | PoseBusters `dock` config on the top pose |
| Rescoring | Crystal pose scored by both: intermolecular energy (the scoring functions alone) and affinity |

Both programs' output PDBQT is read the same way, with Meeko, which rebuilds
each pose's molecule from the `REMARK SMILES` lines the file carries.

## Regression suite

`tests/parity.rs` checks parity on every CI run. It does not need Vina or Python. Its fixtures in
`tests/parity_data/` are 20 complexes from the PoseBusters set. The PoseBusters set comes from the
PDB, so it can be redistributed; PDBbind cannot. The complexes are chosen to cover ligand size,
torsion count, halogens, phosphorus, sulfur, metal sites and Meeko macrocycles. Each receptor is
trimmed to whole residues near what is scored. `golden.tsv` holds Vina 1.2.7's scores for the
trimmed files.

| Test | Checks |
|---|---|
| `crystal_scores_match_vina`, `crystal_scores_match_vinardo` | Crystal-pose intermolecular energy and affinity within 0.02 kcal/mol of Vina's |
| `vina_poses_rescore` | Vina's top docked pose, read from its output, within 0.02 of Vina's own score |
| `docking_finds_crystal_pose` (`--ignored`, release) | Docking 5 small complexes (seed 42): top-1 heavy-atom RMSD < 2 Å and affinity no more than 0.3 above Vina's |

Today rustdock is within 0.016 kcal/mol on every fixture. A 1.7% change to the H-bond weight, or a
7.5 Å cutoff, makes the Vina-scoring tests fail. Docking takes about 8 s on 4 threads.

To regenerate the fixtures after a PoseBusters run (`parity.py run --subset posebusters`):

```bash
benchmark/parity/.venv/bin/python benchmark/parity/make_golden.py
```

The script checks that trimming leaves Vina's scores unchanged. It also checks that each docked
complex's in-order RMSD equals the symmetry-aware RMSD, which is what the Rust test relies on.
