# Docking parity: rustdock-vina vs AutoDock Vina 1.2.7

CASF-2016 core set, subset `pilot` (29 complexes). Identical inputs for both: Meeko-prepared ligand from an RDKit conformer, Open Babel receptor, box = crystal ligand extent + 10 A, exhaustiveness 8, 9 modes, seed 42, 16 CPUs. Success = top pose within 2 A heavy-atom RMSD (symmetry-aware) of the crystal ligand.

| | AutoDock Vina | rustdock-vina |
|---|---|---|
| Docked without error | 29/29 | 29/29 |
| Top-1 success (< 2 A) | 16/29 | 11/29 |
| Any-mode success (< 2 A) | 24/29 | 16/29 |
| Median top-1 RMSD (A) | 1.31 | 3.94 |
| Top-1 PoseBusters-valid | 26/29 | 27/29 |
| Median wall time (s) | 2.4 | 0.5 |

Crystal-pose rescoring, rustdock minus Vina:

- intermolecular energy (scoring function alone): median -0.00, median |diff| 0.00, max |diff| 0.01 kcal/mol
- affinity (with torsional normalization): median -0.00, median |diff| 0.00, max |diff| 0.01 kcal/mol

| Code | Heavy | Rot | Vina RMSD | Rust RMSD | Vina best | Rust best | Vina aff | Rust aff | Xtal inter Vina | Xtal inter Rust | Vina PB | Rust PB | Vina s | Rust s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5c28 | 12 | 3 | 4.74 | 4.71 | 0.98 | 2.61 | -5.94 | -6.12 | -6.01 | -6.01 | yes | yes | 1.2 | 0.2 |
| 4m0y | 22 | 3 | 7.25 | 6.65 | 4.25 | 6.05 | -9.17 | -9.18 | -8.52 | -8.52 | yes | yes | 2.5 | 0.5 |
| 5a7b | 31 | 10 | 1.94 | 15.81 | 1.94 | 7.58 | -7.01 | -4.94 | -10.26 | -10.26 | yes | yes | 10.0 | 0.8 |
| 3k5v | 27 | 5 | 3.34 | 2.07 | 2.14 | 2.07 | -9.82 | -8.61 | -11.43 | -11.43 | yes | yes | 4.4 | 1.0 |
| 3o9i | 42 | 16 | 3.91 | 10.98 | 1.93 | 6.21 | -9.35 | -8.09 | -17.74 | -17.75 | no | yes | 26.5 | 1.2 |
| 3dx1 | 9 | 4 | 2.07 | 2.22 | 2.07 | 1.88 | -5.84 | -5.91 | -5.81 | -5.82 | yes | yes | 1.3 | 0.3 |
| 2vvn | 16 | 6 | 2.96 | 6.28 | 1.14 | 5.79 | -7.30 | -6.74 | -8.50 | -8.49 | yes | yes | 2.0 | 0.7 |
| 4eo8 | 27 | 7 | 1.01 | 4.22 | 0.97 | 4.22 | -7.04 | -5.01 | -10.15 | -10.15 | yes | yes | 5.9 | 0.9 |
| 3f3e | 9 | 6 | 0.78 | 0.52 | 0.57 | 0.52 | -6.67 | -6.80 | -7.67 | -7.67 | yes | yes | 1.1 | 0.3 |
| 3udh | 14 | 0 | 0.37 | 0.39 | 0.37 | 0.39 | -7.51 | -7.54 | -7.11 | -7.12 | yes | yes | 0.9 | 0.2 |
| 1r5y | 13 | 2 | 0.17 | 0.16 | 0.17 | 0.16 | -7.75 | -7.77 | -7.34 | -7.35 | yes | yes | 1.2 | 0.3 |
| 4gkm | 20 | 5 | 4.29 | 5.59 | 0.76 | 5.24 | -8.65 | -8.95 | -8.99 | -8.98 | yes | yes | 2.4 | 0.5 |
| 2hb1 | 14 | 4 | 4.19 | 4.37 | 1.26 | 0.63 | -6.33 | -5.97 | -6.87 | -6.87 | yes | yes | 1.6 | 0.3 |
| 5tmn | 32 | 17 | 2.20 | 8.61 | 2.20 | 5.26 | -8.34 | -5.15 | -12.15 | -12.14 | no | yes | 11.9 | 1.3 |
| 3pxf | 21 | 3 | 0.87 | 0.86 | 0.87 | 0.68 | -10.63 | -10.61 | -11.82 | -11.82 | yes | yes | 2.0 | 0.4 |
| 4jia | 31 | 7 | 1.31 | 1.28 | 0.97 | 1.28 | -9.78 | -9.49 | -11.97 | -11.97 | yes | yes | 4.7 | 0.8 |
| 3acw | 21 | 3 | 1.27 | 1.44 | 1.27 | 1.36 | -10.88 | -10.93 | -10.96 | -10.95 | yes | yes | 2.4 | 0.5 |
| 3u9q | 12 | 9 | 1.16 | 1.44 | 1.16 | 1.25 | -5.96 | -5.62 | -7.56 | -7.56 | yes | yes | 2.2 | 0.3 |
| 1owh | 24 | 5 | 3.17 | 9.57 | 0.87 | 1.55 | -8.62 | -7.72 | -9.30 | -9.30 | yes | yes | 3.7 | 0.7 |
| 3b27 | 16 | 4 | 0.48 | 0.42 | 0.48 | 0.42 | -7.18 | -7.20 | -7.67 | -7.67 | yes | yes | 1.6 | 0.3 |
| 1p1n | 15 | 5 | 0.99 | 1.00 | 0.99 | 0.38 | -8.33 | -8.29 | -8.88 | -8.87 | yes | yes | 1.6 | 0.3 |
| 3fv1 | 21 | 8 | 0.62 | 3.12 | 0.62 | 2.83 | -7.32 | -7.39 | -13.03 | -13.03 | yes | yes | 3.1 | 0.6 |
| 2y5h | 31 | 9 | 2.44 | 7.02 | 1.01 | 5.16 | -8.45 | -6.36 | -10.31 | -10.32 | yes | yes | 6.2 | 0.8 |
| 1uto | 9 | 3 | 1.04 | 3.94 | 1.01 | 0.83 | -4.53 | -4.41 | -4.79 | -4.79 | yes | yes | 1.1 | 0.2 |
| 4j21 | 20 | 3 | 0.24 | 0.23 | 0.24 | 0.23 | -11.84 | -11.87 | -12.69 | -12.69 | yes | yes | 1.9 | 0.3 |
| 3ozs | 33 | 10 | 1.03 | 7.75 | 1.03 | 7.00 | -9.82 | -6.28 | -14.61 | -14.62 | yes | no | 11.3 | 1.3 |
| 3n76 | 22 | 8 | 0.66 | 0.92 | 0.66 | 0.44 | -11.00 | -10.82 | -14.14 | -14.14 | no | no | 3.7 | 1.2 |
| 3nq9 | 10 | 7 | 2.74 | 6.82 | 2.05 | 1.71 | -4.41 | -4.85 | -3.74 | -3.74 | yes | yes | 1.4 | 0.2 |
| 1o0h | 27 | 9 | 2.97 | 11.60 | 1.25 | 8.07 | -6.48 | -6.17 | -7.68 | -7.68 | yes | yes | 6.7 | 0.6 |

## PoseBusters failures (top pose)

- `3o9i` vina: internal_steric_clash,internal_energy
- `5tmn` vina: internal_steric_clash,internal_energy
- `3ozs` rustdock: double_bond_flatness
- `3n76` vina: double_bond_flatness
- `3n76` rustdock: double_bond_flatness
