#!/usr/bin/env python
"""Docking parity between rustdock-vina and AutoDock Vina on CASF-2016.

Both dockers get identical inputs: the same receptor PDBQT, the same ligand
PDBQT (prepared by Meeko from an RDKit conformer, so the crystal pose never
leaks into the search), the same box, exhaustiveness, mode count, seed and
CPU count. Poses are scored against the crystal ligand by symmetry-aware
heavy-atom RMSD and checked with PoseBusters; the crystal pose is also
rescored by both programs to separate scoring differences from search ones.

Usage (from the repository root):
    benchmark/parity/.venv/bin/python benchmark/parity/parity.py run \\
        --pdbbind /path/to/PDBbind --subset pilot --cpu 16
    benchmark/parity/.venv/bin/python benchmark/parity/parity.py report

Every step caches its outputs under --work, so an interrupted run resumes and
`report` can be rerun without docking again.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
VINA = HERE / "bin" / "vina"
RUSTDOCK = REPO / "target" / "release" / "vina"
CORESET = HERE / "casf2016_coreset.tsv"

EXHAUSTIVENESS = 8
NUM_MODES = 9
SEED = 42
BOX_PADDING = 10.0  # Angstroms added to the crystal ligand's extent per axis
SUCCESS_RMSD = 2.0

RDLogger.DisableLog("rdApp.*")


# --------------------------------------------------------------------------
# Selection
# --------------------------------------------------------------------------


def select(subset: str, data: Path) -> list[str]:
    """Complex codes for a subset.

    CASF-2016 core set: full is all 285; targets the best-resolved complex of
    each of the 57 targets; pilot those of the odd-numbered targets (29).
    posebusters: every complex directory under `data` (the PoseBusters
    benchmark set).
    """
    if subset == "posebusters":
        return sorted(p.name for p in data.iterdir() if p.is_dir())
    core = pd.read_csv(CORESET, sep="\t", comment="#", header=None,
                       names=["code", "resolution", "year", "target"])
    if subset == "full":
        return sorted(core.code)
    best = (core.sort_values(["target", "resolution", "code"])
                .groupby("target").head(1))
    if subset == "pilot":
        best = best[best.target % 2 == 1]
    return list(best.sort_values("target").code)


def find_entry(pdbbind: Path, code: str) -> Path:
    """The directory holding `<code>_protein.pdb`: `<root>/<code>` (PoseBusters)
    or `<root>/<group>/<code>` (PDBbind)."""
    for d in [pdbbind / code, *sorted(pdbbind.glob(f"*/{code}"))]:
        if d.parent.name != "demo" and (d / f"{code}_protein.pdb").exists():
            return d
    raise FileNotFoundError(f"{code} not found under {pdbbind}")


# --------------------------------------------------------------------------
# Preparation
# --------------------------------------------------------------------------


def load_crystal_ligand(entry: Path, code: str) -> Chem.Mol:
    """The crystal ligand with explicit hydrogens, from SDF or else MOL2."""
    mol = Chem.MolFromMolFile(str(entry / f"{code}_ligand.sdf"), removeHs=False)
    if mol is None:
        mol = Chem.MolFromMol2File(str(entry / f"{code}_ligand.mol2"), removeHs=False)
    if mol is None:
        raise ValueError("RDKit cannot read the crystal ligand")
    return Chem.AddHs(mol, addCoords=True)


def meeko_pdbqt(mol: Chem.Mol) -> str:
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    setups = MoleculePreparation().prepare(mol)
    text, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise ValueError(f"Meeko could not write PDBQT: {err}")
    return text


def prepare(code: str, entry: Path, out: Path) -> dict:
    """Receptor and ligand PDBQTs plus the search box for one complex."""
    meta_path = out / "prep.json"
    if meta_path.exists():
        return json.loads(meta_path.read_text())
    out.mkdir(parents=True, exist_ok=True)

    # Receptor: protein without waters, typed and charged by Open Babel.
    # PDBbind proteins carry hydrogens; PoseBusters' do not, and without
    # polar hydrogens a receptor has no H-bond donors, so add them at pH 7.4.
    protein = out / "protein_dry.pdb"
    lines = (entry / f"{code}_protein.pdb").read_text().splitlines()
    protein.write_text("\n".join(l for l in lines if l[17:20] != "HOH") + "\n")
    atoms = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
    has_h = any(l[76:78].strip() == "H" for l in atoms)
    protonate = [] if has_h else ["-p", "7.4"]
    subprocess.run(["obabel", "-ipdb", str(protein), *protonate, "-opdbqt", "-xr",
                    "-O", str(out / "receptor.pdbqt")],
                   check=True, capture_output=True)

    # Ligand: the crystal pose (for rescoring) and an RDKit conformer of the
    # same molecule (for docking), both with Meeko, so atom order matches.
    xtal = load_crystal_ligand(entry, code)
    Chem.MolToMolFile(xtal, str(out / "ligand_xtal.sdf"))
    (out / "ligand_xtal.pdbqt").write_text(meeko_pdbqt(xtal))

    start = Chem.Mol(xtal)
    start.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    if AllChem.EmbedMolecule(start, params) != 0:
        raise ValueError("RDKit could not embed a starting conformer")
    AllChem.MMFFOptimizeMolecule(start)
    (out / "ligand_start.pdbqt").write_text(meeko_pdbqt(start))

    heavy = np.array([xtal.GetConformer().GetAtomPosition(a.GetIdx())
                      for a in xtal.GetAtoms() if a.GetAtomicNum() > 1])
    center = heavy.mean(axis=0)
    size = heavy.max(axis=0) - heavy.min(axis=0) + BOX_PADDING
    meta = {
        "code": code,
        "center": [round(float(c), 3) for c in center],
        "size": [round(float(s), 3) for s in size],
        "heavy_atoms": int(len(heavy)),
        "rotatable_bonds": int(Chem.rdMolDescriptors.CalcNumRotatableBonds(xtal)),
    }
    meta_path.write_text(json.dumps(meta, indent=1))
    return meta


# --------------------------------------------------------------------------
# Docking and rescoring
# --------------------------------------------------------------------------


def dock(tool: str, out: Path, meta: dict, cpu: int, scoring: str) -> dict:
    result_path = out / f"{tool}_dock.json"
    if result_path.exists():
        return json.loads(result_path.read_text())
    poses = out / f"{tool}_out.pdbqt"
    c, s = meta["center"], meta["size"]
    rec, lig = str(out / "receptor.pdbqt"), str(out / "ligand_start.pdbqt")
    env = None
    if tool == "vina":
        cmd = [str(VINA), "--receptor", rec, "--ligand", lig,
               "--center_x", str(c[0]), "--center_y", str(c[1]), "--center_z", str(c[2]),
               "--size_x", str(s[0]), "--size_y", str(s[1]), "--size_z", str(s[2]),
               "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", str(NUM_MODES),
               "--seed", str(SEED), "--cpu", str(cpu), "--scoring", scoring,
               "--out", str(poses)]
    else:
        cmd = [str(RUSTDOCK), "dock", "--receptor", rec, "--ligand", lig,
               # `=` keeps a negative first coordinate from reading as a flag.
               "--center=" + ",".join(map(str, c)), "--size=" + ",".join(map(str, s)),
               "--exhaustiveness", str(EXHAUSTIVENESS), "--num-modes", str(NUM_MODES),
               "--seed", str(SEED), "--scoring", scoring, "--out", str(poses)]
        env = {"RAYON_NUM_THREADS": str(cpu), "RUST_LOG": "warn",
               "PATH": "/usr/bin:/bin"}
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    seconds = time.perf_counter() - t0
    result = {"ok": proc.returncode == 0 and poses.exists(), "seconds": round(seconds, 2)}
    if not result["ok"]:
        result["error"] = (proc.stderr or proc.stdout).strip()[-2000:]
    result_path.write_text(json.dumps(result, indent=1))
    return result


def rescore(tool: str, out: Path, scoring: str) -> dict:
    """Affinity and intermolecular energy the tool assigns to the crystal pose.

    Intermolecular energy compares the scoring functions alone; affinity also
    includes each program's torsional normalization.
    """
    rec, lig = str(out / "receptor.pdbqt"), str(out / "ligand_xtal.pdbqt")
    if tool == "vina":
        cmd = [str(VINA), "--receptor", rec, "--ligand", lig, "--score_only",
               "--autobox", "--cpu", "1", "--scoring", scoring]
        patterns = {"affinity": r"Estimated Free Energy of Binding\s*:\s*(-?[\d.]+)",
                    "inter": r"\(1\) Final Intermolecular Energy\s*:\s*(-?[\d.]+)"}
    else:
        cmd = [str(RUSTDOCK), "score", "--receptor", rec, "--ligand", lig,
               "--scoring", scoring]
        patterns = {"affinity": r"Affinity:\s*(-?[\d.]+)",
                    "inter": r"Intermolecular energy:\s*(-?[\d.]+)"}
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          env={"RUST_LOG": "warn", "PATH": "/usr/bin:/bin"})
    found = {}
    for key, pattern in patterns.items():
        m = re.search(pattern, proc.stdout)
        found[f"{tool}_xtal_{key}"] = float(m.group(1)) if m else None
    return found


# --------------------------------------------------------------------------
# Pose evaluation
# --------------------------------------------------------------------------

ATOM_RE = re.compile(r"^(ATOM  |HETATM)")


def model_coords(pdbqt: str) -> list[tuple[float | None, np.ndarray]]:
    """(reported affinity, atom coordinates) for each MODEL of a docking output."""
    models, coords, energy = [], [], None
    for line in pdbqt.splitlines():
        if line.startswith("MODEL"):
            coords, energy = [], None
        elif line.startswith("REMARK VINA RESULT:"):
            energy = float(line.split()[3])
        elif ATOM_RE.match(line):
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        elif line.startswith("ENDMDL"):
            models.append((energy, np.array(coords)))
    return models


def pose_mols(tool: str, out: Path) -> list[tuple[float | None, Chem.Mol]]:
    """Docked poses as RDKit molecules with bond orders, best first.

    Both programs' output is read the same way, with Meeko: it keeps the
    ligand's REMARK SMILES and torsion tree, so the file alone rebuilds the
    molecule.
    """
    from meeko import PDBQTMolecule, RDKitMolCreate

    text = (out / f"{tool}_out.pdbqt").read_text()
    models = model_coords(text)
    mol = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(text, skip_typing=True))[0]
    poses = []
    for i, (energy, _) in enumerate(models):
        single = Chem.Mol(mol)
        conf = Chem.Conformer(mol.GetConformers()[i])
        single.RemoveAllConformers()
        single.AddConformer(conf, assignId=True)
        poses.append((energy, single))
    return poses


def rmsd(pose: Chem.Mol, xtal: Chem.Mol) -> float:
    """Symmetry-aware heavy-atom RMSD, in place (no superposition)."""
    return rdMolAlign.CalcRMS(Chem.RemoveHs(pose), Chem.RemoveHs(xtal))


def posebusters_valid(pose: Chem.Mol, out: Path) -> tuple[bool | None, str]:
    """Whether the pose passes every PoseBusters docking check, and which fail."""
    from posebusters import PoseBusters

    pred = out / "_pb_pose.sdf"
    Chem.MolToMolFile(pose, str(pred))
    try:
        frame = PoseBusters(config="dock").bust(
            mol_pred=pred, mol_true=out / "ligand_xtal.sdf",
            mol_cond=out / "protein_dry.pdb", full_report=False)
    except Exception as exc:  # noqa: BLE001 - report, do not abort the table
        return None, f"error: {exc}"[:80]
    checks = frame.iloc[0]
    checks = checks[[c for c in checks.index if isinstance(checks[c], (bool, np.bool_))]]
    failed = [c for c, ok in checks.items() if not ok]
    return not failed, ",".join(failed)


def vina_scores_top_pose(out: Path, scoring: str) -> float | None:
    """Vina's affinity for rustdock's top pose, read from rustdock's output.

    Checks both that Vina can read the file and that rustdock reports the
    affinity Vina would.
    """
    lines = (out / "rustdock_out.pdbqt").read_text().splitlines()
    model = []
    for line in lines[1:]:
        if line.startswith("ENDMDL"):
            break
        model.append(line)
    top = out / "_rustdock_top.pdbqt"
    top.write_text("\n".join(model) + "\n")
    proc = subprocess.run([str(VINA), "--receptor", str(out / "receptor.pdbqt"),
                           "--ligand", str(top), "--score_only", "--autobox", "--cpu", "1",
                           "--scoring", scoring],
                          capture_output=True, text=True)
    m = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?[\d.]+)", proc.stdout)
    return float(m.group(1)) if m else None


def evaluate(tool: str, out: Path, docked: dict, scoring: str) -> dict:
    row = {f"{tool}_seconds": docked["seconds"]}
    if not docked["ok"]:
        row[f"{tool}_error"] = docked.get("error", "")[:200]
        return row
    xtal = Chem.MolFromMolFile(str(out / "ligand_xtal.sdf"), removeHs=False)
    try:
        poses = pose_mols(tool, out)
    except Exception as exc:  # noqa: BLE001
        row[f"{tool}_error"] = f"pose export: {exc}"[:200]
        return row
    rmsds = [rmsd(mol, xtal) for _, mol in poses]
    valid, failed = posebusters_valid(poses[0][1], out)
    row.update({
        f"{tool}_affinity": poses[0][0],
        f"{tool}_rmsd_top1": round(rmsds[0], 3),
        f"{tool}_rmsd_best": round(min(rmsds), 3),
        f"{tool}_modes": len(poses),
        f"{tool}_pb_valid": valid,
        f"{tool}_pb_failed": failed,
    })
    if tool == "rustdock":
        row["rustdock_top_vina_score"] = vina_scores_top_pose(out, scoring)
    return row


# --------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------


def run(args: argparse.Namespace) -> None:
    if not VINA.exists():
        sys.exit(f"missing {VINA}; see benchmark/parity/README.md")
    subprocess.run(["cargo", "build", "--release", "--quiet"], cwd=REPO, check=True)
    codes = args.codes or select(args.subset, args.pdbbind)
    rows = []
    for n, code in enumerate(codes, 1):
        out = args.work / code
        row = {"code": code}
        try:
            meta = prepare(code, find_entry(args.pdbbind, code), out)
            row.update(heavy_atoms=meta["heavy_atoms"], rotatable_bonds=meta["rotatable_bonds"])
            for tool in ("vina", "rustdock"):
                docked = dock(tool, out, meta, args.cpu, args.scoring)
                row.update(evaluate(tool, out, docked, args.scoring))
                row.update(rescore(tool, out, args.scoring))
        except Exception as exc:  # noqa: BLE001 - one bad entry must not stop the run
            row["prep_error"] = str(exc)[:200]
        rows.append(row)
        print(f"[{n}/{len(codes)}] {code}: "
              f"vina {row.get('vina_rmsd_top1', '-')} A {row.get('vina_seconds', '-')} s | "
              f"rustdock {row.get('rustdock_rmsd_top1', '-')} A {row.get('rustdock_seconds', '-')} s",
              flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(args.work / "results.csv", index=False)
    write_report(frame, args)


def summarize(frame: pd.DataFrame, tool: str) -> dict:
    top1 = frame.get(f"{tool}_rmsd_top1")
    best = frame.get(f"{tool}_rmsd_best")
    docked = top1.notna() if top1 is not None else pd.Series(False, index=frame.index)
    n = len(frame)
    return {
        "docked": f"{int(docked.sum())}/{n}",
        "top1": f"{int((top1 < SUCCESS_RMSD).sum())}/{n}" if top1 is not None else "-",
        "best": f"{int((best < SUCCESS_RMSD).sum())}/{n}" if best is not None else "-",
        "median_rmsd": f"{top1.median():.2f}" if top1 is not None else "-",
        "pb_valid": (f"{int(frame[f'{tool}_pb_valid'].fillna(False).astype(bool).sum())}/{n}"
                     if f"{tool}_pb_valid" in frame else "-"),
        "median_seconds": f"{frame[f'{tool}_seconds'].median():.1f}"
        if f"{tool}_seconds" in frame else "-",
    }


def top_pose_check(frame: pd.DataFrame) -> str:
    if "rustdock_top_vina_score" not in frame:
        return "no complexes (not run)"
    read = frame["rustdock_top_vina_score"].notna()
    d = (frame["rustdock_affinity"] - frame["rustdock_top_vina_score"]).abs()[read]
    text = f"{int(read.sum())}/{len(frame)} complexes"
    if read.any():
        text += (f"; its score differs from rustdock's reported affinity by median "
                 f"{d.median():.3f}, max {d.max():.3f} kcal/mol")
    return text


def fmt(value, digits=2) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "-"
    if isinstance(value, (bool, np.bool_)):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def write_report(frame: pd.DataFrame, args: argparse.Namespace) -> None:
    vina, rust = summarize(frame, "vina"), summarize(frame, "rustdock")
    def gap(key: str) -> str:
        cols = {f"rustdock_xtal_{key}", f"vina_xtal_{key}"}
        if not cols <= set(frame):
            return "unavailable"
        d = frame[f"rustdock_xtal_{key}"] - frame[f"vina_xtal_{key}"]
        if not d.notna().any():
            return "unavailable"
        return (f"median {d.median():+.2f}, median |diff| {d.abs().median():.2f}, "
                f"max |diff| {d.abs().max():.2f} kcal/mol")
    lines = [
        f"# Docking parity: rustdock-vina vs AutoDock Vina 1.2.7 ({args.scoring} scoring)",
        "",
        (f"PoseBusters benchmark set ({len(frame)} complexes). "
         if args.subset == "posebusters" else
         f"CASF-2016 core set, subset `{args.subset}` ({len(frame)} complexes). ")
        +
        f"Identical inputs for both: Meeko-prepared ligand from an RDKit conformer, "
        f"Open Babel receptor, box = crystal ligand extent + {BOX_PADDING:g} A, "
        f"exhaustiveness {EXHAUSTIVENESS}, {NUM_MODES} modes, seed {SEED}, "
        f"{args.cpu} CPUs. Success = top pose within {SUCCESS_RMSD:g} A heavy-atom "
        "RMSD (symmetry-aware) of the crystal ligand.",
        "",
        "| | AutoDock Vina | rustdock-vina |",
        "|---|---|---|",
        f"| Docked without error | {vina['docked']} | {rust['docked']} |",
        f"| Top-1 success (< 2 A) | {vina['top1']} | {rust['top1']} |",
        f"| Any-mode success (< 2 A) | {vina['best']} | {rust['best']} |",
        f"| Median top-1 RMSD (A) | {vina['median_rmsd']} | {rust['median_rmsd']} |",
        f"| Top-1 PoseBusters-valid | {vina['pb_valid']} | {rust['pb_valid']} |",
        f"| Median wall time (s) | {vina['median_seconds']} | {rust['median_seconds']} |",
        "",
        "Crystal-pose rescoring, rustdock minus Vina:",
        "",
        f"- intermolecular energy (scoring function alone): {gap('inter')}",
        f"- affinity (with torsional normalization): {gap('affinity')}",
        "",
        "Output round trip: Vina read rustdock's top pose from its output file in "
        + top_pose_check(frame) + ".",
        "",
        "| Code | Heavy | Rot | Vina RMSD | Rust RMSD | Vina best | Rust best | "
        "Vina aff | Rust aff | Xtal inter Vina | Xtal inter Rust | Vina PB | Rust PB | "
        "Vina s | Rust s |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for _, r in frame.iterrows():
        g = r.get
        lines.append("| " + " | ".join([
            str(r["code"]), fmt(g("heavy_atoms")), fmt(g("rotatable_bonds")),
            fmt(g("vina_rmsd_top1")), fmt(g("rustdock_rmsd_top1")),
            fmt(g("vina_rmsd_best")), fmt(g("rustdock_rmsd_best")),
            fmt(g("vina_affinity")), fmt(g("rustdock_affinity")),
            fmt(g("vina_xtal_inter")), fmt(g("rustdock_xtal_inter")),
            fmt(g("vina_pb_valid")), fmt(g("rustdock_pb_valid")),
            fmt(g("vina_seconds"), 1), fmt(g("rustdock_seconds"), 1),
        ]) + " |")
    errors = [(r["code"], c, r[c]) for _, r in frame.iterrows()
              for c in frame.columns if c.endswith("error") and isinstance(r[c], str)]
    if errors:
        lines += ["", "## Errors", ""]
        lines += [f"- `{code}` {col}: {msg}" for code, col, msg in errors]
    failures = [(r["code"], t, r[f"{t}_pb_failed"]) for _, r in frame.iterrows()
                for t in ("vina", "rustdock")
                if isinstance(r.get(f"{t}_pb_failed"), str) and r[f"{t}_pb_failed"]]
    if failures:
        lines += ["", "## PoseBusters failures (top pose)", ""]
        lines += [f"- `{code}` {tool}: {checks}" for code, tool, checks in failures]
    suffix = "" if args.scoring == "vina" else f"_{args.scoring}"
    report = args.report or (HERE / f"PARITY_{args.subset}{suffix}.md")
    report.write_text("\n".join(lines) + "\n")
    print(f"wrote {report}")


def report(args: argparse.Namespace) -> None:
    write_report(pd.read_csv(args.work / "results.csv"), args)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name, func in (("run", run), ("report", report)):
        p = sub.add_parser(name)
        p.set_defaults(func=func)
        p.add_argument("--work", type=Path,
                       help="work directory (default: work, or work-<scoring>)")
        p.add_argument("--scoring", choices=["vina", "vinardo"], default="vina")
        p.add_argument("--subset", choices=["pilot", "targets", "full", "posebusters"],
                       default="pilot")
        p.add_argument("--cpu", type=int, default=16)
        p.add_argument("--report", type=Path)
    sub.choices["run"].add_argument("--pdbbind", "--data", type=Path, required=True,
                                    help="dataset root: PDBbind (<group>/<code>/) or "
                                         "the PoseBusters benchmark set (<code>/)")
    sub.choices["run"].add_argument("--codes", nargs="*", help="override the subset")
    args = parser.parse_args()
    if args.work is None:
        name = "work" if args.subset != "posebusters" else "work-posebusters"
        if args.scoring != "vina":
            name += f"-{args.scoring}"
        args.work = HERE / name
    args.work.mkdir(parents=True, exist_ok=True)
    args.func(args)


if __name__ == "__main__":
    main()
