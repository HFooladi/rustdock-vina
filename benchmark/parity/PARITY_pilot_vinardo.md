# Docking parity: rustdock-vina vs AutoDock Vina 1.2.7 (vinardo scoring)

CASF-2016 core set, subset `pilot` (29 complexes). Identical inputs for both: Meeko-prepared ligand from an RDKit conformer, Open Babel receptor, box = crystal ligand extent + 10 A, exhaustiveness 8, 9 modes, seed 42, 16 CPUs. Success = top pose within 2 A heavy-atom RMSD (symmetry-aware) of the crystal ligand.

| | AutoDock Vina | rustdock-vina |
|---|---|---|
| Docked without error | 29/29 | 29/29 |
| Top-1 success (< 2 A) | 18/29 | 19/29 |
| Any-mode success (< 2 A) | 21/29 | 25/29 |
| Median top-1 RMSD (A) | 1.06 | 1.08 |
| Top-1 PoseBusters-valid | 26/29 | 27/29 |
| Median wall time (s) | 2.4 | 1.7 |

Crystal-pose rescoring, rustdock minus Vina:

- intermolecular energy (scoring function alone): median +0.00, median |diff| 0.00, max |diff| 0.02 kcal/mol
- affinity (with torsional normalization): median +0.00, median |diff| 0.00, max |diff| 0.01 kcal/mol

Output round trip: Vina read rustdock's top pose from its output file in 29/29 complexes; its score differs from rustdock's reported affinity by median 0.003, max 0.007 kcal/mol.

| Code | Heavy | Rot | Vina RMSD | Rust RMSD | Vina best | Rust best | Vina aff | Rust aff | Xtal inter Vina | Xtal inter Rust | Vina PB | Rust PB | Vina s | Rust s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5c28 | 12 | 3 | 4.78 | 2.53 | 0.49 | 1.81 | -4.57 | -4.90 | -4.69 | -4.69 | yes | yes | 1.2 | 0.5 |
| 4m0y | 22 | 3 | 7.32 | 7.32 | 5.19 | 5.60 | -6.72 | -6.72 | -5.65 | -5.65 | yes | yes | 2.6 | 1.9 |
| 5a7b | 31 | 10 | 1.53 | 1.53 | 1.53 | 1.53 | -5.76 | -5.75 | -8.26 | -8.25 | yes | yes | 12.0 | 11.3 |
| 3k5v | 27 | 5 | 3.34 | 3.34 | 2.10 | 2.28 | -7.31 | -7.31 | -8.35 | -8.35 | yes | yes | 5.2 | 4.2 |
| 3o9i | 42 | 16 | 0.80 | 0.78 | 0.80 | 0.78 | -7.17 | -7.16 | -13.96 | -13.96 | yes | yes | 36.8 | 28.5 |
| 3dx1 | 9 | 4 | 2.33 | 2.25 | 2.15 | 2.25 | -5.08 | -5.05 | -6.07 | -6.07 | yes | yes | 1.3 | 0.7 |
| 2vvn | 16 | 6 | 6.28 | 6.28 | 2.54 | 1.61 | -5.06 | -5.06 | -6.41 | -6.41 | yes | yes | 2.4 | 1.6 |
| 4eo8 | 27 | 7 | 0.79 | 0.80 | 0.79 | 0.80 | -5.85 | -5.85 | -8.59 | -8.59 | yes | yes | 6.2 | 4.5 |
| 3f3e | 9 | 6 | 0.48 | 0.48 | 0.48 | 0.48 | -6.90 | -6.90 | -7.81 | -7.81 | yes | yes | 1.1 | 0.5 |
| 3udh | 14 | 0 | 0.42 | 0.39 | 0.42 | 0.39 | -6.44 | -6.45 | -5.86 | -5.86 | yes | yes | 0.9 | 0.5 |
| 1r5y | 13 | 2 | 0.22 | 0.22 | 0.22 | 0.22 | -6.57 | -6.57 | -5.74 | -5.74 | yes | yes | 1.2 | 0.7 |
| 4gkm | 20 | 5 | 5.59 | 5.59 | 0.48 | 0.52 | -6.72 | -6.72 | -6.45 | -6.45 | yes | yes | 2.5 | 1.8 |
| 2hb1 | 14 | 4 | 1.45 | 1.37 | 1.45 | 1.37 | -4.38 | -4.40 | -4.66 | -4.66 | yes | yes | 1.9 | 1.3 |
| 5tmn | 32 | 17 | 2.15 | 2.76 | 1.88 | 1.82 | -7.25 | -7.64 | -9.63 | -9.63 | yes | yes | 14.0 | 11.8 |
| 3pxf | 21 | 3 | 0.77 | 0.77 | 0.77 | 0.77 | -9.49 | -9.52 | -9.50 | -9.50 | yes | yes | 2.1 | 1.4 |
| 4jia | 31 | 7 | 1.06 | 1.08 | 1.04 | 1.08 | -7.20 | -7.21 | -8.24 | -8.23 | yes | yes | 5.0 | 5.6 |
| 3acw | 21 | 3 | 1.05 | 1.07 | 1.05 | 1.07 | -8.46 | -8.47 | -8.88 | -8.88 | yes | yes | 2.2 | 1.3 |
| 3u9q | 12 | 9 | 0.81 | 0.75 | 0.81 | 0.75 | -5.93 | -5.83 | -7.37 | -7.37 | yes | yes | 2.1 | 1.2 |
| 1owh | 24 | 5 | 2.18 | 2.18 | 2.14 | 0.89 | -6.84 | -6.83 | -6.66 | -6.66 | yes | yes | 4.6 | 2.9 |
| 3b27 | 16 | 4 | 0.53 | 0.55 | 0.53 | 0.55 | -5.02 | -5.02 | -5.10 | -5.10 | yes | yes | 1.6 | 1.0 |
| 1p1n | 15 | 5 | 0.88 | 0.88 | 0.88 | 0.88 | -6.63 | -6.63 | -6.10 | -6.10 | yes | yes | 1.7 | 0.9 |
| 3fv1 | 21 | 8 | 0.62 | 0.62 | 0.62 | 0.62 | -5.65 | -5.65 | -10.27 | -10.26 | no | yes | 3.4 | 2.6 |
| 2y5h | 31 | 9 | 2.31 | 2.32 | 2.31 | 2.32 | -5.70 | -5.67 | -6.10 | -6.11 | yes | yes | 7.1 | 6.4 |
| 1uto | 9 | 3 | 1.04 | 1.02 | 1.04 | 1.02 | -3.90 | -3.89 | -3.50 | -3.50 | yes | yes | 1.1 | 0.6 |
| 4j21 | 20 | 3 | 0.38 | 0.39 | 0.38 | 0.39 | -9.73 | -9.75 | -9.29 | -9.28 | yes | yes | 2.1 | 1.7 |
| 3ozs | 33 | 10 | 6.17 | 6.15 | 3.61 | 0.89 | -6.49 | -6.45 | -11.03 | -11.05 | no | no | 15.8 | 12.7 |
| 3n76 | 22 | 8 | 0.84 | 0.84 | 0.84 | 0.84 | -9.70 | -9.71 | -12.33 | -12.32 | no | no | 3.7 | 3.1 |
| 3nq9 | 10 | 7 | 1.67 | 1.91 | 1.67 | 1.35 | -3.79 | -3.91 | -3.23 | -3.23 | yes | yes | 1.2 | 0.7 |
| 1o0h | 27 | 9 | 2.01 | 1.88 | 2.01 | 1.88 | -4.68 | -4.81 | -5.75 | -5.75 | yes | yes | 7.4 | 6.2 |

## PoseBusters failures (top pose)

- `3fv1` vina: internal_energy
- `3ozs` vina: double_bond_flatness
- `3ozs` rustdock: double_bond_flatness
- `3n76` vina: double_bond_flatness
- `3n76` rustdock: double_bond_flatness
