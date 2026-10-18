# Docking parity: rustdock-vina vs AutoDock Vina 1.2.7

CASF-2016 core set, subset `pilot` (29 complexes). Identical inputs for both: Meeko-prepared ligand from an RDKit conformer, Open Babel receptor, box = crystal ligand extent + 10 A, exhaustiveness 8, 9 modes, seed 42, 16 CPUs. Success = top pose within 2 A heavy-atom RMSD (symmetry-aware) of the crystal ligand.

| | AutoDock Vina | rustdock-vina |
|---|---|---|
| Docked without error | 29/29 | 29/29 |
| Top-1 success (< 2 A) | 16/29 | 8/29 |
| Any-mode success (< 2 A) | 24/29 | 12/29 |
| Median top-1 RMSD (A) | 1.31 | 4.41 |
| Top-1 PoseBusters-valid | 26/29 | 27/29 |
| Median wall time (s) | 2.4 | 0.5 |

Crystal-pose rescoring, rustdock minus Vina:

- intermolecular energy (scoring function alone): median +1.73, median |diff| 1.73, max |diff| 6.57 kcal/mol
- affinity (with torsional normalization): median +1.29, median |diff| 1.29, max |diff| 5.07 kcal/mol

| Code | Heavy | Rot | Vina RMSD | Rust RMSD | Vina best | Rust best | Vina aff | Rust aff | Xtal inter Vina | Xtal inter Rust | Vina PB | Rust PB | Vina s | Rust s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5c28 | 12 | 3 | 4.74 | 2.53 | 0.98 | 2.53 | -5.94 | -5.80 | -6.01 | -6.01 | yes | yes | 1.2 | 0.2 |
| 4m0y | 22 | 3 | 7.25 | 6.21 | 4.25 | 5.87 | -9.17 | -8.29 | -8.52 | -6.47 | yes | yes | 2.5 | 0.5 |
| 5a7b | 31 | 10 | 1.94 | 13.01 | 1.94 | 9.03 | -7.01 | -4.94 | -10.26 | -10.26 | yes | yes | 10.0 | 0.8 |
| 3k5v | 27 | 5 | 3.34 | 3.36 | 2.14 | 2.02 | -9.82 | -9.53 | -11.43 | -11.43 | yes | yes | 4.4 | 1.0 |
| 3o9i | 42 | 16 | 3.91 | 5.23 | 1.93 | 4.22 | -9.35 | -7.69 | -17.74 | -16.01 | no | yes | 26.5 | 1.3 |
| 3dx1 | 9 | 4 | 2.07 | 3.43 | 2.07 | 2.05 | -5.84 | -4.71 | -5.81 | -3.33 | yes | yes | 1.3 | 0.3 |
| 2vvn | 16 | 6 | 2.96 | 3.00 | 1.14 | 3.00 | -7.30 | -6.30 | -8.50 | -7.09 | yes | yes | 2.0 | 0.7 |
| 4eo8 | 27 | 7 | 1.01 | 8.21 | 0.97 | 3.69 | -7.04 | -5.51 | -10.15 | -8.67 | yes | yes | 5.9 | 1.0 |
| 3f3e | 9 | 6 | 0.78 | 1.96 | 0.57 | 0.98 | -6.67 | -5.44 | -7.67 | -3.84 | yes | yes | 1.1 | 0.3 |
| 3udh | 14 | 0 | 0.37 | 4.70 | 0.37 | 0.42 | -7.51 | -6.71 | -7.11 | -5.08 | yes | yes | 0.9 | 0.2 |
| 1r5y | 13 | 2 | 0.17 | 5.06 | 0.17 | 3.06 | -7.75 | -5.65 | -7.34 | -3.69 | yes | yes | 1.2 | 0.3 |
| 4gkm | 20 | 5 | 4.29 | 6.41 | 0.76 | 6.08 | -8.65 | -7.35 | -8.99 | -6.32 | yes | yes | 2.4 | 0.5 |
| 2hb1 | 14 | 4 | 4.19 | 4.41 | 1.26 | 3.67 | -6.33 | -5.30 | -6.87 | -4.72 | yes | yes | 1.6 | 0.3 |
| 5tmn | 32 | 17 | 2.20 | 6.37 | 2.20 | 6.01 | -8.34 | -5.32 | -12.15 | -9.96 | no | yes | 11.9 | 1.3 |
| 3pxf | 21 | 3 | 0.87 | 1.13 | 0.87 | 0.99 | -10.63 | -9.79 | -11.82 | -10.30 | yes | yes | 2.0 | 0.4 |
| 4jia | 31 | 7 | 1.31 | 10.28 | 0.97 | 1.77 | -9.78 | -8.50 | -11.97 | -11.01 | yes | yes | 4.7 | 0.9 |
| 3acw | 21 | 3 | 1.27 | 1.42 | 1.27 | 1.36 | -10.88 | -10.66 | -10.96 | -10.95 | yes | yes | 2.4 | 0.6 |
| 3u9q | 12 | 9 | 1.16 | 1.40 | 1.16 | 1.40 | -5.96 | -5.42 | -7.56 | -6.44 | yes | yes | 2.2 | 0.3 |
| 1owh | 24 | 5 | 3.17 | 10.15 | 0.87 | 5.06 | -8.62 | -5.93 | -9.30 | -6.30 | yes | yes | 3.7 | 0.8 |
| 3b27 | 16 | 4 | 0.48 | 0.60 | 0.48 | 0.60 | -7.18 | -6.87 | -7.67 | -7.35 | yes | yes | 1.6 | 0.3 |
| 1p1n | 15 | 5 | 0.99 | 3.74 | 0.99 | 1.57 | -8.33 | -6.82 | -8.88 | -5.28 | yes | yes | 1.6 | 0.3 |
| 3fv1 | 21 | 8 | 0.62 | 3.15 | 0.62 | 2.98 | -7.32 | -5.96 | -13.03 | -6.45 | yes | yes | 3.1 | 0.6 |
| 2y5h | 31 | 9 | 2.44 | 9.25 | 1.01 | 5.98 | -8.45 | -6.12 | -10.31 | -10.11 | yes | yes | 6.2 | 0.8 |
| 1uto | 9 | 3 | 1.04 | 1.51 | 1.01 | 1.23 | -4.53 | -3.75 | -4.79 | -3.99 | yes | yes | 1.1 | 0.2 |
| 4j21 | 20 | 3 | 0.24 | 0.46 | 0.24 | 0.46 | -11.84 | -11.11 | -12.69 | -11.76 | yes | yes | 1.9 | 0.4 |
| 3ozs | 33 | 10 | 1.03 | 15.15 | 1.03 | 9.39 | -9.82 | -5.92 | -14.61 | -12.47 | yes | no | 11.3 | 1.3 |
| 3n76 | 22 | 8 | 0.66 | 0.69 | 0.66 | 0.69 | -11.00 | -8.48 | -14.14 | -10.65 | no | no | 3.7 | 1.2 |
| 3nq9 | 10 | 7 | 2.74 | 6.82 | 2.05 | 1.48 | -4.41 | -4.85 | -3.74 | -3.74 | yes | yes | 1.4 | 0.2 |
| 1o0h | 27 | 9 | 2.97 | 8.64 | 1.25 | 6.51 | -6.48 | -4.90 | -7.68 | -4.19 | yes | yes | 6.7 | 0.6 |

## PoseBusters failures (top pose)

- `3o9i` vina: internal_steric_clash,internal_energy
- `5tmn` vina: internal_steric_clash,internal_energy
- `3ozs` rustdock: double_bond_flatness
- `3n76` vina: double_bond_flatness
- `3n76` rustdock: double_bond_flatness
