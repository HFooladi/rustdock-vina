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
