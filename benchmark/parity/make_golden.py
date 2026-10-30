#!/usr/bin/env python
"""Build the Vina parity fixtures that `tests/parity.rs` checks in CI.

Takes a stratified sample of about 20 complexes from a finished PoseBusters
parity run (`parity.py run --subset posebusters`). It trims each receptor to the
residues near the search box, and records AutoDock Vina 1.2.7's scores for them
in `tests/parity_data/golden.tsv`:

- the crystal pose's intermolecular energy and affinity, with Vina and Vinardo
  scoring;
- the same for Vina's top docked pose (Vina scoring), with that pose stored as
  `vina_top.pdbqt`;
- for a few small complexes that both programs dock within 2 A, Vina's docked
  affinity and top-1 RMSD, so the test can check that rustdock's search still
  finds the crystal pose;
- for about ten complexes, the crystal pose scored with flexible side chains:
  up to three residues next to the ligand, prepared with Meeko as
  `flex_rigid.pdbqt` and `flex.pdbqt`; those also docked are docked with the
  side chains too, recording Vina's affinity and top-1 RMSD.

Only PoseBusters complexes are used. They come from the PDB, so they can be
redistributed, which PDBbind/CASF structures cannot.

Usage (from the repository root, after a PoseBusters parity run):
    benchmark/parity/.venv/bin/python benchmark/parity/make_golden.py
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
VINA = HERE / "bin" / "vina"

# Receptor residues are kept whole when any atom is within the 8 A cutoff of
# a scored ligand atom, plus this much slack.
TRIM_SLACK = 1.0

METALS = {"Zn", "Mg", "Fe", "Mn", "Ca", "Co", "Ni", "Cu", "Na", "K"}
ATOM_RECORD = ("ATOM", "HETATM")

# (label, predicate on a feature row). Each one picks the first unused complex,
# in code order, that matches.
STRATA = [
    ("tiny", lambda r: r.ha <= 12),
    ("small rigid", lambda r: 13 <= r.ha <= 20 and r.rot <= 2),
    ("medium", lambda r: 21 <= r.ha <= 30 and 3 <= r.rot <= 6),
    ("large", lambda r: 31 <= r.ha <= 40 and 7 <= r.rot <= 12),
    ("very flexible", lambda r: r.ha > 45 and r.rot > 15),
    ("no torsions", lambda r: r.rot == 0),
    ("fluorine", lambda r: "F" in r.types),
    ("chlorine", lambda r: "Cl" in r.types),
    ("phosphorus", lambda r: "P" in r.types),
    ("sulfur", lambda r: "S" in r.types or "SA" in r.types),
    ("zinc site", lambda r: "Zn" in r.metals),
    ("magnesium site", lambda r: "Mg" in r.metals),
    ("iron site", lambda r: "Fe" in r.metals),
    ("macrocycle", lambda r: r.macro),
    ("macrocycle", lambda r: r.macro),
]
# Docked in CI, so small enough to run in seconds, and found by both programs
# with room to spare.
DOCK_COUNT = 5

MK_PREPARE_RECEPTOR = HERE / ".venv" / "bin" / "mk_prepare_receptor.py"
# Complexes scored with flexible side chains, and how many residues each.
FLEX_COUNT = 10
FLEX_RESIDUES = 3
# A residue is flexible when a side-chain heavy atom is this close to the
# crystal ligand.
FLEX_CONTACT = 4.0
# Side chains with a torsion; Meeko knows their templates.
FLEX_KINDS = {"ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE", "LEU", "LYS",
              "MET", "PHE", "SER", "THR", "TRP", "TYR", "VAL"}
BACKBONE = {"N", "CA", "C", "O", "OXT"}


def dock_candidate(r, work: Path) -> bool:
    if not (r.ha <= 20 and r.rot <= 4 and r.vina_rmsd < 1.5 and r.rustdock_rmsd < 1.5):
        return False
    # The Rust test measures RMSD in atom order; that is only right when it
    # agrees with the symmetry-aware RMSD of the parity run.
    d = work / r.code
    top = first_model((d / "vina_out.pdbqt").read_text())
    return abs(heavy_rmsd((d / "ligand_xtal.pdbqt").read_text(), top) - r.vina_rmsd) < 0.01


def atoms(pdbqt: str) -> list[str]:
    return [l for l in pdbqt.splitlines() if l.startswith(ATOM_RECORD)]


def xyz(line: str) -> np.ndarray:
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def ad_type(line: str) -> str:
    return line[77:79].strip()


def first_model(pdbqt: str) -> str:
    """The first MODEL of a docking output, without the MODEL/ENDMDL lines."""
    lines = pdbqt.splitlines()
    end = next(i for i, l in enumerate(lines) if l.startswith("ENDMDL"))
    return "\n".join(l for l in lines[:end] if not l.startswith("MODEL")) + "\n"


def heavy_rmsd(a: str, b: str) -> float:
    """In-order heavy-atom RMSD of two PDBQTs of the same ligand."""
    ha = [l for l in atoms(a) if ad_type(l) not in ("H", "HD")]
    hb = [l for l in atoms(b) if ad_type(l) not in ("H", "HD")]
    if [ad_type(l) for l in ha] != [ad_type(l) for l in hb]:
        return float("nan")
    d = np.array([xyz(l) for l in ha]) - np.array([xyz(l) for l in hb])
    return float(np.sqrt((d**2).sum(axis=1).mean()))


def features(work: Path) -> pd.DataFrame:
    frame = pd.read_csv(work / "results.csv")
    rows = []
    for r in frame.itertuples():
        if isinstance(r.vina_error, str) or pd.isna(r.vina_xtal_inter):
            continue
        d = work / r.code
        types = {ad_type(l) for l in atoms((d / "ligand_xtal.pdbqt").read_text())}
        rec = {ad_type(l) for l in atoms((d / "receptor.pdbqt").read_text())}
        rows.append(dict(code=r.code, ha=r.heavy_atoms, rot=r.rotatable_bonds,
                         types=types, metals=rec & METALS,
                         macro=bool(types & {"G0", "CG0"}),
                         vina_rmsd=r.vina_rmsd_top1, rustdock_rmsd=r.rustdock_rmsd_top1,
                         vina_affinity=r.vina_affinity))
    return pd.DataFrame(rows).sort_values("code").reset_index(drop=True)


def select(feats: pd.DataFrame, work: Path) -> tuple[list[tuple[str, str]], list[str]]:
    """(code, stratum) for the scoring set, and the codes docked in CI."""
    chosen: dict[str, str] = {}
    for label, pred in STRATA:
        for r in feats.itertuples():
            if r.code not in chosen and pred(r):
                chosen[r.code] = label
                break
        else:
            print(f"  no complex for stratum {label!r}")
    dock = []
    for r in feats.itertuples():
        if len(dock) == DOCK_COUNT:
            break
        if dock_candidate(r, work):
            dock.append(r.code)
            chosen.setdefault(r.code, "docked")
    return list(chosen.items()), dock


def trim_receptor(receptor: str, keep_atom) -> str:
    """Whole residues with at least one atom for which `keep_atom(xyz)` holds.

    Keeping whole residues keeps every bond, and with them Vina's atom typing.
    """
    keep = {l[17:27] for l in atoms(receptor) if keep_atom(xyz(l))}  # residue id
    return "".join(l + "\n" for l in receptor.splitlines()
                   if not l.startswith(ATOM_RECORD) or l[17:27] in keep)


def near_ligands(*ligands: str):
    """Within the cutoff of an atom of any of these ligand PDBQTs."""
    points = np.array([xyz(l) for lig in ligands for l in atoms(lig)])
    reach = 8.0 + TRIM_SLACK
    return lambda p: bool((((points - p) ** 2).sum(axis=1) < reach**2).any())


def near_box(center: np.ndarray, size: np.ndarray, ligand: str):
    """Everything docking can reach: the box plus the ligand's radius and the
    cutoff, as rustdock's own receptor pruning keeps."""
    coords = np.array([xyz(l) for l in atoms(ligand)])
    radius = np.sqrt(((coords - coords.mean(axis=0)) ** 2).sum(axis=1)).max()
    half = size / 2 + radius + 8.0 + TRIM_SLACK
    return lambda p: bool(np.all(np.abs(p - center) <= half))


def contact_residues(protein: Path, ligand: str) -> list[str]:
    """Up to FLEX_RESIDUES residues ("A:79") whose side chains touch the ligand,
    closest first."""
    points = np.array([xyz(l) for l in atoms(ligand) if ad_type(l) not in ("H", "HD")])
    found = []
    for chain in gemmi.read_structure(str(protein))[0]:
        for res in chain:
            if res.name not in FLEX_KINDS:
                continue
            side = [a for a in res if a.name not in BACKBONE and a.element.name != "H"]
            if not side:
                continue
            d = min(np.sqrt(((points - np.array(a.pos.tolist())) ** 2).sum(axis=1)).min()
                    for a in side)
            if d < FLEX_CONTACT:
                found.append((d, f"{chain.name}:{res.seqid.num}"))
    return [r for _, r in sorted(found)[:FLEX_RESIDUES]]


def prepare_flex(src: Path, dst: Path, code: str, ligand: str) -> list[str]:
    """Meeko's rigid and flexible receptor files for the residues next to the
    ligand, as `flex_rigid.pdbqt` and `flex.pdbqt` in `dst`; [] on failure."""
    # The protein without the ligand's own residue, which Meeko would try to
    # template.
    ccd = code.split("_")[1]
    lines = (src / "protein_dry.pdb").read_text().splitlines()
    protein = dst / "protein_for_flex.pdb"
    protein.write_text("\n".join(l for l in lines if l[17:20].strip() != ccd) + "\n")
    residues = contact_residues(protein, ligand)
    if len(residues) < 2:
        protein.unlink()
        return []
    base = dst / "flex_prep"
    proc = subprocess.run([str(MK_PREPARE_RECEPTOR), "--read_pdb", str(protein), "-o", str(base),
                           "-p", "-a", "-f", ",".join(residues)], capture_output=True, text=True)
    protein.unlink()
    rigid, flex = Path(f"{base}_rigid.pdbqt"), Path(f"{base}_flex.pdbqt")
    if proc.returncode != 0 or not flex.exists():
        print(f"  {code}: Meeko could not make flex residues {residues}")
        for f in (rigid, flex):
            f.unlink(missing_ok=True)
        return []
    rigid.rename(dst / "flex_rigid.pdbqt")
    flex.rename(dst / "flex.pdbqt")
    return residues


def flex_box(*pdbqts: str) -> tuple[np.ndarray, np.ndarray]:
    """A box around every atom of these PDBQTs with 4 A to spare: Vina scores
    only when the ligand and every flexible atom are inside its grid."""
    points = np.array([xyz(l) for text in pdbqts for l in atoms(text)])
    lo, hi = points.min(axis=0), points.max(axis=0)
    return np.round((lo + hi) / 2, 3), np.round(hi - lo + 8.0, 3)


def vina_score_flex(rigid: Path, flex: Path, ligand: Path, scoring: str,
                    center: np.ndarray, size: np.ndarray) -> dict[str, float]:
    """`vina --score_only --flex` energy terms.

    Vina's affinity subtracts an unbound energy whose flex-receptor term comes
    from the grid, while its internal energy uses exact pairs, so the two do
    not cancel. The second run sets `--unbound_energy` to the internal energy,
    which makes the affinity the normalized intermolecular energy alone.
    """
    def run(*extra: str) -> str:
        box = [f"--{k}_{a}" for k in ("center", "size") for a in "xyz"]
        values = [str(v) for v in (*center, *size)]
        return subprocess.run([str(VINA), "--receptor", str(rigid), "--flex", str(flex),
                               "--ligand", str(ligand), "--score_only", "--cpu", "1",
                               "--scoring", scoring,
                               *[x for pair in zip(box, values) for x in pair], *extra],
                              capture_output=True, text=True, check=True).stdout

    def term(text: str, label: str) -> float:
        return float(re.search(re.escape(label) + r"\s*:\s*(-?[\d.]+)", text).group(1))

    out = run()
    intra = term(out, "(2) Final Total Internal Energy")
    # Vina treats an unbound energy of exactly 0 as unset.
    unbound = intra if intra != 0 else 1e-9
    consistent = run("--unbound_energy", str(unbound))
    return {
        "inter": term(out, "(1) Final Intermolecular Energy"),
        "affinity": term(consistent, "Estimated Free Energy of Binding"),
        "lig_flex": term(out, "Ligand - Flex side chains"),
        "flex_rec": term(out, "Flex   - Receptor"),
        "flex_flex": term(out, "Flex   - Flex side chains"),
    }


def add_flex(row: dict, src: Path, dst: Path) -> bool:
    """Prepare flexible residues for one complex and add Vina's flex scores of
    its crystal pose to `row`, and for a docked complex Vina's flex docking;
    False if it has too few contact residues."""
    xtal = (src / "ligand_xtal.pdbqt").read_text()
    residues = prepare_flex(src, dst, row["code"], xtal)
    if not residues:
        return False
    flex_text = (dst / "flex.pdbqt").read_text()
    full = dst / "flex_rigid_full.pdbqt"
    (dst / "flex_rigid.pdbqt").rename(full)
    center, size = flex_box(xtal, flex_text)
    if row["dock"]:
        keep = near_box(center, size, (src / "ligand_start.pdbqt").read_text())
    else:
        keep = near_ligands(xtal, flex_text)
    (dst / "flex_rigid.pdbqt").write_text(trim_receptor(full.read_text(), keep))
    row["flex_residues"] = " ".join(residues)
    for axis, c, sz in zip("xyz", center, size):
        row[f"flex_center_{axis}"], row[f"flex_size_{axis}"] = c, sz
    for scoring in ("vina", "vinardo"):
        rest = (dst / "flex.pdbqt", dst / "ligand_xtal.pdbqt", scoring, center, size)
        terms = vina_score_flex(dst / "flex_rigid.pdbqt", *rest)
        # Trimming must not change any score: check against the full receptor.
        full_terms = vina_score_flex(full, *rest)
        for name in ("inter", "flex_rec"):
            assert abs(terms[name] - full_terms[name]) < 1e-3, (row["code"], scoring, name)
        for name, value in terms.items():
            row[f"{scoring}_flex_{name}"] = value
    full.unlink()
    if row["dock"]:
        affinity, rmsd = vina_dock_flex(dst, center, size)
        row["vina_flex_dock_affinity"], row["vina_flex_dock_rmsd"] = affinity, rmsd
        # Docked in CI only if Vina finds the crystal pose with room to spare.
        row["flex_dock"] = rmsd < 1.5
    return True


def vina_dock_flex(dst: Path, center: np.ndarray, size: np.ndarray) -> tuple[float, float]:
    """Vina's top affinity and in-order top-1 heavy-atom RMSD to the crystal
    pose, docking `ligand_start.pdbqt` with the flexible side chains (seed
    42, exhaustiveness 8)."""
    out = dst / "vina_flex_out.pdbqt"
    box = [x for k, values in (("center", center), ("size", size))
           for a, v in zip("xyz", values) for x in (f"--{k}_{a}", str(v))]
    subprocess.run([str(VINA), "--receptor", str(dst / "flex_rigid.pdbqt"),
                    "--flex", str(dst / "flex.pdbqt"), "--ligand", str(dst / "ligand_start.pdbqt"),
                    "--seed", "42", "--exhaustiveness", "8", "--cpu", "16", *box,
                    "--out", str(out)], capture_output=True, text=True, check=True)
    top = first_model(out.read_text())
    out.unlink()
    affinity = float(re.search(r"REMARK VINA RESULT:\s*(-?[\d.]+)", top).group(1))
    ligand = top[:top.index("BEGIN_RES")]
    return affinity, heavy_rmsd((dst / "ligand_xtal.pdbqt").read_text(), ligand)


def vina_score(receptor: Path, ligand: Path, scoring: str) -> tuple[float, float]:
    """(intermolecular energy, affinity) from `vina --score_only`."""
    proc = subprocess.run([str(VINA), "--receptor", str(receptor), "--ligand", str(ligand),
                           "--score_only", "--autobox", "--cpu", "1", "--scoring", scoring],
                          capture_output=True, text=True, check=True)
    inter = re.search(r"\(1\) Final Intermolecular Energy\s*:\s*(-?[\d.]+)", proc.stdout)
    aff = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?[\d.]+)", proc.stdout)
    return float(inter.group(1)), float(aff.group(1))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--work", type=Path, default=HERE / "work-posebusters",
                        help="work directory of a finished PoseBusters parity run")
    parser.add_argument("--out", type=Path, default=REPO / "tests" / "parity_data")
    args = parser.parse_args()

    version = subprocess.run([str(VINA), "--version"], capture_output=True,
                             text=True).stdout.strip()
    print(f"Vina: {version}")
    feats = features(args.work)
    picks, dock = select(feats, args.work)
    by_code = feats.set_index("code")

    if args.out.exists():
        shutil.rmtree(args.out)
    rows = []
    for code, label in picks:
        src, dst = args.work / code, args.out / code
        dst.mkdir(parents=True)
        meta = pd.read_json(src / "prep.json", typ="series")
        center, size = np.array(meta["center"]), np.array(meta["size"])

        for name in ("ligand_xtal.pdbqt", "ligand_start.pdbqt"):
            shutil.copy(src / name, dst / name)
        vina_top = first_model((src / "vina_out.pdbqt").read_text())
        (dst / "vina_top.pdbqt").write_text(vina_top)
        if code in dock:
            keep = near_box(center, size, (src / "ligand_start.pdbqt").read_text())
        else:
            keep = near_ligands((src / "ligand_xtal.pdbqt").read_text(), vina_top)
        full = (src / "receptor.pdbqt").read_text()
        (dst / "receptor.pdbqt").write_text(trim_receptor(full, keep))

        row = dict(code=code, stratum=label,
                   center_x=center[0], center_y=center[1], center_z=center[2],
                   size_x=size[0], size_y=size[1], size_z=size[2])
        for scoring in ("vina", "vinardo"):
            inter, aff = vina_score(dst / "receptor.pdbqt", dst / "ligand_xtal.pdbqt", scoring)
            # Trimming must not change any score: check against the full receptor.
            full_inter, _ = vina_score(src / "receptor.pdbqt", src / "ligand_xtal.pdbqt", scoring)
            assert abs(inter - full_inter) < 1e-3, (code, scoring, inter, full_inter)
            row[f"{scoring}_xtal_inter"], row[f"{scoring}_xtal_affinity"] = inter, aff
        row["vina_top_inter"], row["vina_top_affinity"] = vina_score(
            dst / "receptor.pdbqt", dst / "vina_top.pdbqt", "vina")

        row["dock"] = code in dock
        if row["dock"]:
            row["vina_dock_affinity"] = by_code.loc[code, "vina_affinity"]
            row["vina_dock_rmsd"] = by_code.loc[code, "vina_rmsd"]

        rows.append(row)

    # Flexible side chains, docked complexes first.
    flex_done = 0
    for row in sorted(rows, key=lambda r: not r["dock"]):
        if flex_done == FLEX_COUNT:
            break
        if add_flex(row, args.work / row["code"], args.out / row["code"]):
            flex_done += 1

    golden = pd.DataFrame(rows)
    golden.to_csv(args.out / "golden.tsv", sep="\t", index=False, float_format="%.3f")
    (args.out / "README.md").write_text(
        "Generated by `benchmark/parity/make_golden.py` from the PoseBusters "
        f"benchmark set (PDB-derived) with {version}. Do not edit by hand; "
        "see benchmark/parity/README.md.\n")
    total = sum(p.stat().st_size for p in args.out.rglob("*") if p.is_file())
    print(golden[["code", "stratum", "dock", "vina_xtal_affinity",
                  "vinardo_xtal_affinity", "flex_residues"]].to_string(index=False))
    print(f"{len(golden)} complexes, {golden.dock.sum()} docked, "
          f"{golden.flex_residues.notna().sum()} with flex, {total / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
