"""Virtual-screening parity: does rustdock rank actives above decoys as well
as AutoDock Vina 1.2.7 does?

Targets come from DUDE-Z (Stein et al., J. Chem. Inf. Model. 2021, the
successor of DUD-E): a prepared receptor, its crystal ligand, actives
(ChEMBL, several protomers each) and 50 property-matched decoys per active
(ZINC). Each target's actives are docked with a random sample of its decoys,
by both programs with the same inputs, box and settings, single-threaded on
pinned cores. A compound scores its best affinity over its protomers.

Metrics per target: ROC AUC, BEDROC (alpha 20) and enrichment factors at
1% and 5%. rustdock's difference from Vina is set against Vina's own
seed-to-seed difference, the noise floor of any comparison of two runs.

Usage (from the repository root, with the parity venv):
  P=benchmark/parity/.venv/bin/python
  $P -I benchmark/screening/screen.py run --work DIR [--targets ...] [--jobs 64]
  $P -I benchmark/screening/screen.py report --work DIR --report FILE
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import random
import subprocess
import sys
import time
import urllib.request
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
VINA = REPO / "benchmark" / "parity" / "bin" / "vina"
RUSTDOCK = REPO / "target" / "release" / "rustdock"
DUDEZ = "https://dudez.docking.org"

# Ten diverse DUDE-Z targets: two GPCR/nuclear-receptor pockets each, two
# kinases, and five enzymes of different classes.
TARGETS = ["AA2AR", "ABL1", "ACES", "AMPC", "ANDR", "EGFR", "HIVPR", "HS90A",
           "PARP1", "THRB"]
DECOYS = None  # decoys sampled per target; None docks them all
SEED = 42  # decoy sample, RDKit conformers and the main docking seed
SEEDS = (42, 1)  # docking seeds: the second gives each program's noise floor
EXHAUSTIVENESS = 8
NUM_MODES = 9
BOX_PADDING = 10.0  # Angstroms added to the crystal ligand's extent per axis
TOOLS = ("vina", "rustdock")
BOOTSTRAP = 2_000


# --- preparation -----------------------------------------------------------

def fetch(url: str, path: Path) -> Path:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(url, timeout=120) as response:
            data = response.read()
        path.write_bytes(data)
    return path


def read_smi(path: Path) -> list[tuple[str, str]]:
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            rows.append((parts[0], parts[1]))
    return rows


def ligand_pdbqt(smiles: str) -> str:
    """A Meeko PDBQT of an RDKit conformer (ETKDG + MMFF) of `smiles`."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit cannot parse the SMILES")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError("RDKit could not embed a conformer")
    AllChem.MMFFOptimizeMolecule(mol)
    setups = MoleculePreparation().prepare(mol)
    text, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise ValueError(f"Meeko could not write PDBQT: {err}")
    return text


def prepare_target(target: str, work: Path, decoys_n: int | None) -> dict:
    """Receptor PDBQT, box and ligand list (actives and sampled decoys)."""
    out = work / target
    meta_path = out / "target.json"
    if meta_path.exists():
        return json.loads(meta_path.read_text())
    grids = f"{DUDEZ}/DOCKING_GRIDS_AND_POSES/{target}"
    sets = f"{DUDEZ}/DUDE-Z-benchmark-grids/{target}"
    rec = fetch(f"{grids}/rec.crg.pdb", out / "src" / "rec.crg.pdb")
    xtal = fetch(f"{grids}/xtal-lig.pdb", out / "src" / "xtal-lig.pdb")
    actives = read_smi(fetch(f"{sets}/ligands.smi", out / "src" / "ligands.smi"))
    decoys = read_smi(fetch(f"{sets}/decoys.smi", out / "src" / "decoys.smi"))

    # DUDE-Z's receptor carries its polar hydrogens; Open Babel types it.
    subprocess.run(["obabel", "-ipdb", str(rec), "-opdbqt", "-xr",
                    "-O", str(out / "receptor.pdbqt")], check=True, capture_output=True)

    heavy = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                      for l in xtal.read_text().splitlines()
                      if l.startswith(("ATOM", "HETATM")) and l[76:78].strip() != "H"])
    center = heavy.mean(axis=0)
    size = heavy.max(axis=0) - heavy.min(axis=0) + BOX_PADDING

    ids = sorted({cid for _, cid in decoys})
    chosen = set(ids if decoys_n is None else random.Random(SEED).sample(ids, min(decoys_n, len(ids))))
    ligands = []
    for label, rows in ((1, actives), (0, [r for r in decoys if r[1] in chosen])):
        counts: dict[str, int] = {}
        for smiles, cid in rows:
            k = counts.get(cid, 0)
            counts[cid] = k + 1
            ligands.append({"name": f"{cid}_{k}", "compound": cid, "active": label,
                            "smiles": smiles})
    meta = {"target": target,
            "center": [round(float(c), 3) for c in center],
            "size": [round(float(s), 3) for s in size],
            "ligands": ligands}
    meta_path.write_text(json.dumps(meta, indent=1))
    return meta


def prepare_ligand(job: tuple[str, str, str]) -> tuple[str, str, str | None]:
    out, name, smiles = job
    path = Path(out) / "ligands" / f"{name}.pdbqt"
    error_path = path.with_suffix(".error")
    if path.exists():
        return out, name, None
    if error_path.exists():
        return out, name, error_path.read_text()
    try:
        text = ligand_pdbqt(smiles)
    except Exception as exc:  # noqa: BLE001 - any failure leaves the ligand out
        error_path.write_text(str(exc))
        return out, name, str(exc)
    path.write_text(text)
    return out, name, None


# --- docking ---------------------------------------------------------------

def pin(slots) -> None:
    os.sched_setaffinity(0, [slots.get()])


def top_affinity(pdbqt: Path) -> float | None:
    for line in pdbqt.read_text().splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
    return None


def dock(job: tuple) -> dict:
    """Dock one ligand with one tool and seed, single-threaded."""
    out, name, tool, seed, center, size = job
    out = Path(out)
    result_path = out / "docked" / f"{name}.{tool}.s{seed}.json"
    if result_path.exists():
        return json.loads(result_path.read_text())
    poses = out / "docked" / f"{name}.{tool}.s{seed}.pdbqt"
    rec, lig = str(out / "receptor.pdbqt"), str(out / "ligands" / f"{name}.pdbqt")
    c, s = center, size
    if tool == "vina":
        cmd = [str(VINA), "--receptor", rec, "--ligand", lig,
               "--center_x", str(c[0]), "--center_y", str(c[1]), "--center_z", str(c[2]),
               "--size_x", str(s[0]), "--size_y", str(s[1]), "--size_z", str(s[2]),
               "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", str(NUM_MODES),
               "--seed", str(seed), "--cpu", "1", "--out", str(poses)]
        env = None
    else:
        cmd = [str(RUSTDOCK), "dock", "--receptor", rec, "--ligand", lig,
               "--center=" + ",".join(map(str, c)), "--size=" + ",".join(map(str, s)),
               "--exhaustiveness", str(EXHAUSTIVENESS), "--num-modes", str(NUM_MODES),
               "--seed", str(seed), "--cpu", "1", "--out", str(poses)]
        env = {"RAYON_NUM_THREADS": "1", "RUST_LOG": "warn", "PATH": "/usr/bin:/bin"}
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    seconds = time.perf_counter() - t0
    affinity = top_affinity(poses) if proc.returncode == 0 and poses.exists() else None
    result = {"name": name, "tool": tool, "seed": seed, "affinity": affinity,
              "seconds": round(seconds, 2)}
    if affinity is None:
        result["error"] = (proc.stderr or proc.stdout).strip()[-2000:]
    result_path.write_text(json.dumps(result))
    return result


def run(args: argparse.Namespace) -> None:
    if not VINA.exists():
        sys.exit(f"missing {VINA}; see benchmark/parity/README.md")
    subprocess.run(["cargo", "build", "--release", "--quiet"], cwd=REPO, check=True)
    metas = [prepare_target(t, args.work, args.decoys) for t in args.targets]

    prep = []
    for meta in metas:
        out = args.work / meta["target"]
        (out / "ligands").mkdir(exist_ok=True)
        (out / "docked").mkdir(exist_ok=True)
        prep += [(str(out), l["name"], l["smiles"]) for l in meta["ligands"]
                 if args.limit is None or i_within(l, meta, args.limit)]
    failed = set()
    with ProcessPoolExecutor(args.jobs) as pool:
        for out, name, error in pool.map(prepare_ligand, prep, chunksize=8):
            if error:
                failed.add((out, name))
    print(f"prepared {len(prep) - len(failed)} ligands, {len(failed)} failed", flush=True)

    jobs = [(out, name, tool, seed, meta["center"], meta["size"])
            for meta in metas
            for out, name, _ in prep if out == str(args.work / meta["target"])
            if (out, name) not in failed
            for tool in TOOLS for seed in args.seeds]
    # Interleave tools, targets and seeds so both programs see the same load.
    random.Random(SEED).shuffle(jobs)
    slots = mp.Queue()
    for core in range(args.jobs):
        slots.put(core)
    t0 = time.time()
    with ProcessPoolExecutor(args.jobs, initializer=pin, initargs=(slots,)) as pool:
        futures = [pool.submit(dock, job) for job in jobs]
        for n, future in enumerate(as_completed(futures), 1):
            future.result()
            if n % 500 == 0 or n == len(jobs):
                rate = n / max(time.time() - t0, 1e-9)
                print(f"[{n}/{len(jobs)}] docked, {(len(jobs) - n) / rate / 3600:.1f} h left",
                      flush=True)
    report(args)


def i_within(ligand: dict, meta: dict, limit: int) -> bool:
    """For pilots: the first `limit` actives and decoys of a target."""
    same = [l for l in meta["ligands"] if l["active"] == ligand["active"]]
    return same.index(ligand) < limit


# --- metrics ---------------------------------------------------------------

def roc_auc(scores: np.ndarray, labels: np.ndarray) -> float:
    """Probability that an active outscores a decoy (lower affinity is better)."""
    from scipy.stats import rankdata

    ranks = rankdata(-scores)  # best (most negative) affinity gets the top rank
    n_act, n_dec = labels.sum(), len(labels) - labels.sum()
    return float((ranks[labels == 1].sum() - n_act * (n_act + 1) / 2) / (n_act * n_dec))


def order(scores: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Labels sorted from best to worst affinity; ties broken pessimistically
    (decoys first), so a tie never helps either program."""
    return labels[np.lexsort((labels, scores))]


def enrichment(scores: np.ndarray, labels: np.ndarray, fraction: float) -> float:
    ranked = order(scores, labels)
    n = max(1, int(round(fraction * len(ranked))))
    return float(ranked[:n].sum() / n / labels.mean())


def bedroc(scores: np.ndarray, labels: np.ndarray, alpha: float = 20.0) -> float:
    """BEDROC of Truchon and Bayly (2007)."""
    ranked = order(scores, labels)
    n, n_act = len(ranked), int(ranked.sum())
    ra = n_act / n
    ri = np.nonzero(ranked)[0] + 1
    s = np.exp(-alpha * ri / n).sum()
    rie = s / (ra * (1 - np.exp(-alpha)) / (np.exp(alpha / n) - 1))
    return float(rie * ra * np.sinh(alpha / 2) / (np.cosh(alpha / 2) - np.cosh(alpha / 2 - alpha * ra))
                 + 1 / (1 - np.exp(alpha * (1 - ra))))


METRICS = {"AUC": roc_auc,
           "BEDROC": bedroc,
           "EF1%": lambda s, l: enrichment(s, l, 0.01),
           "EF5%": lambda s, l: enrichment(s, l, 0.05)}


# --- report ----------------------------------------------------------------

def collect(work: Path, targets: list[str]) -> pd.DataFrame:
    """One row per compound and run: its best affinity over its protomers.
    Only compounds every run docked are kept, so all runs rank one set."""
    rows = []
    for target in targets:
        meta = json.loads((work / target / "target.json").read_text())
        info = {l["name"]: l for l in meta["ligands"]}
        for path in (work / target / "docked").glob("*.json"):
            r = json.loads(path.read_text())
            rows.append({"target": target, "compound": info[r["name"]]["compound"],
                         "active": info[r["name"]]["active"], "run": f"{r['tool']}.s{r['seed']}",
                         "affinity": r["affinity"], "seconds": r["seconds"]})
    frame = pd.DataFrame(rows)
    seconds = frame.groupby("run")["seconds"].mean()
    failures = frame[frame["affinity"].isna()].groupby("run").size()
    best = (frame.dropna(subset=["affinity"])
            .groupby(["target", "compound", "active", "run"], as_index=False)["affinity"].min())
    wide = best.pivot_table(index=["target", "compound", "active"], columns="run",
                            values="affinity").dropna().reset_index()
    wide.attrs["seconds"] = seconds.to_dict()
    wide.attrs["failures"] = failures.to_dict()
    return wide


def metric_table(wide: pd.DataFrame, runs: list[str]) -> pd.DataFrame:
    rows = []
    for target, group in wide.groupby("target"):
        labels = group["active"].to_numpy()
        for run in runs:
            scores = group[run].to_numpy()
            row = {"target": target, "run": run, "actives": int(labels.sum()),
                   "decoys": int(len(labels) - labels.sum())}
            row.update({m: f(scores, labels) for m, f in METRICS.items()})
            rows.append(row)
    return pd.DataFrame(rows)


def paired_bootstrap(wide: pd.DataFrame, a: str, b: str, metric: str) -> tuple[float, float, float]:
    """Mean over targets of metric(a) - metric(b), with a 95% CI from
    resampling compounds within each target."""
    rng = np.random.default_rng(SEED)
    f = METRICS[metric]
    groups = [(g["active"].to_numpy(), g[a].to_numpy(), g[b].to_numpy())
              for _, g in wide.groupby("target")]
    point = np.mean([f(x, l) - f(y, l) for l, x, y in groups])
    draws = np.empty(BOOTSTRAP)
    for k in range(BOOTSTRAP):
        diffs = []
        for l, x, y in groups:
            act, dec = np.nonzero(l)[0], np.nonzero(l == 0)[0]
            idx = np.concatenate([rng.choice(act, len(act)), rng.choice(dec, len(dec))])
            diffs.append(f(x[idx], l[idx]) - f(y[idx], l[idx]))
        draws[k] = np.mean(diffs)
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return float(point), float(lo), float(hi)


def report(args: argparse.Namespace) -> None:
    from scipy.stats import spearmanr

    wide = collect(args.work, args.targets)
    runs = [f"{tool}.s{seed}" for tool in TOOLS for seed in args.seeds
            if f"{tool}.s{seed}" in wide.columns]
    table = metric_table(wide, runs)
    table.to_csv(args.work / "metrics.csv", index=False)
    wide.to_csv(args.work / "scores.csv", index=False)
    main_v, main_r = f"vina.s{args.seeds[0]}", f"rustdock.s{args.seeds[0]}"
    noise = len(args.seeds) > 1
    alt_v = f"vina.s{args.seeds[1]}" if noise else None
    alt_r = f"rustdock.s{args.seeds[1]}" if noise else None

    lines = ["# Virtual-screening parity: rustdock vs AutoDock Vina 1.2.7", ""]
    n_act = int(wide["active"].sum())
    lines += [f"DUDE-Z, {wide['target'].nunique()} targets, {n_act} actives and "
              f"{len(wide) - n_act} decoys ({'all' if args.decoys is None else f'a random {args.decoys}'} per target), docked by both "
              f"programs with the same Meeko PDBQTs, box (crystal ligand + {BOX_PADDING:g} A) "
              f"and settings: Vina scoring, exhaustiveness {EXHAUSTIVENESS}, one thread "
              f"per dock on its own core, seeds {', '.join(map(str, args.seeds))}. "
              "A compound scores the best top-pose affinity over its protomers; "
              "compounds that failed in any run are left out of all.", ""]

    lines += ["## Per target (seed %d)" % args.seeds[0], "",
              "| Target | Actives | Decoys | AUC Vina | AUC rustdock | BEDROC Vina | "
              "BEDROC rustdock | EF1% Vina | EF1% rustdock | Spearman ρ |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    for target, group in wide.groupby("target"):
        v = table[(table.target == target) & (table.run == main_v)].iloc[0]
        r = table[(table.target == target) & (table.run == main_r)].iloc[0]
        rho = spearmanr(group[main_v], group[main_r]).statistic
        lines.append(f"| {target} | {v.actives} | {v.decoys} | {v.AUC:.3f} | {r.AUC:.3f} | "
                     f"{v.BEDROC:.3f} | {r.BEDROC:.3f} | {v['EF1%']:.1f} | {r['EF1%']:.1f} | "
                     f"{rho:.3f} |")
    mean = table.groupby("run")[list(METRICS)].mean()
    lines += [f"| **Mean** | | | {mean.loc[main_v, 'AUC']:.3f} | {mean.loc[main_r, 'AUC']:.3f} | "
              f"{mean.loc[main_v, 'BEDROC']:.3f} | {mean.loc[main_r, 'BEDROC']:.3f} | "
              f"{mean.loc[main_v, 'EF1%']:.1f} | {mean.loc[main_r, 'EF1%']:.1f} | |", ""]

    lines += ["## rustdock − Vina, against Vina's seed-to-seed noise", "",
              "Mean over targets of the paired difference, 95% CI from resampling "
              "compounds within each target. A rustdock − Vina difference inside the "
              "range of Vina − Vina (same program, another seed) is noise.", "",
              "| Metric | rustdock − Vina | Vina − Vina (seeds) | rustdock − rustdock (seeds) |",
              "|---|---|---|---|"]
    for metric in METRICS:
        cells = []
        pairs = [(main_r, main_v)]
        if noise:
            pairs += [(alt_v, main_v), (alt_r, main_r)]
        for a, b in pairs:
            p, lo, hi = paired_bootstrap(wide, a, b, metric)
            cells.append(f"{p:+.3f} ({lo:+.3f} to {hi:+.3f})")
        cells += ["-"] * (3 - len(cells))
        lines.append(f"| {metric} | " + " | ".join(cells) + " |")
    lines.append("")

    lines += ["## Per-compound scores", "",
              "| Pair | Spearman ρ | Mean diff (kcal/mol) | Mean abs diff | Within 0.5 kcal/mol |",
              "|---|---|---|---|---|"]
    pairs = [(main_r, main_v)] + ([(alt_v, main_v), (alt_r, main_r)] if noise else [])
    for a, b in pairs:
        d = wide[a] - wide[b]
        lines.append(f"| {a} − {b} | {spearmanr(wide[a], wide[b]).statistic:.3f} | "
                     f"{d.mean():+.3f} | {d.abs().mean():.3f} | {(d.abs() <= 0.5).mean():.1%} |")
    lines.append("")

    seconds, failures = wide.attrs["seconds"], wide.attrs["failures"]
    lines += ["## Time and failures", "", "| Run | Mean seconds per dock (1 thread) | Failed docks |",
              "|---|---|---|"]
    for run in runs:
        lines.append(f"| {run} | {seconds.get(run, float('nan')):.1f} | {failures.get(run, 0)} |")
    lines.append("")
    args.report.write_text("\n".join(lines))
    print("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["run", "report"])
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--targets", nargs="+", default=TARGETS)
    parser.add_argument("--seeds", nargs="+", type=int, default=list(SEEDS))
    parser.add_argument("--jobs", type=int, default=64, help="parallel docks, one core each")
    parser.add_argument("--decoys", type=int, default=DECOYS,
                        help="decoys sampled per target (default: all)")
    parser.add_argument("--limit", type=int, help="pilot: first N actives and decoys per target")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    args.work = args.work.resolve()
    args.report = (args.report or args.work / "SCREENING.md").resolve()
    args.work.mkdir(parents=True, exist_ok=True)
    run(args) if args.command == "run" else report(args)


if __name__ == "__main__":
    main()
