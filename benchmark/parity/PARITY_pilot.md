# Docking parity: rustdock-vina vs AutoDock Vina 1.2.7 (vina scoring)

CASF-2016 core set, subset `pilot` (29 complexes). Identical inputs for both: Meeko-prepared ligand from an RDKit conformer, Open Babel receptor, box = crystal ligand extent + 10 A, exhaustiveness 8, 9 modes, seed 42, 16 CPUs. Success = top pose within 2 A heavy-atom RMSD (symmetry-aware) of the crystal ligand.

| | AutoDock Vina | rustdock-vina |
|---|---|---|
| Docked without error | 29/29 | 29/29 |
| Top-1 success (< 2 A) | 16/29 | 15/29 |
| Any-mode success (< 2 A) | 24/29 | 26/29 |
| Median top-1 RMSD (A) | 1.31 | 1.34 |
| Top-1 PoseBusters-valid | 26/29 | 26/29 |
| Median wall time (s) | 2.4 | 1.6 |

Crystal-pose rescoring, rustdock minus Vina:

- intermolecular energy (scoring function alone): median -0.00, median |diff| 0.00, max |diff| 0.01 kcal/mol
- affinity (with torsional normalization): median -0.00, median |diff| 0.00, max |diff| 0.01 kcal/mol

Output round trip: Vina read rustdock's top pose from its output file in 29/29 complexes; its score differs from rustdock's reported affinity by median 0.003, max 0.007 kcal/mol.

| Code | Heavy | Rot | Vina RMSD | Rust RMSD | Vina best | Rust best | Vina aff | Rust aff | Xtal inter Vina | Xtal inter Rust | Vina PB | Rust PB | Vina s | Rust s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5c28 | 12 | 3 | 4.74 | 4.73 | 0.98 | 0.43 | -5.94 | -6.10 | -6.01 | -6.01 | yes | yes | 1.2 | 0.6 |
| 4m0y | 22 | 3 | 7.25 | 7.28 | 4.25 | 5.24 | -9.17 | -9.16 | -8.52 | -8.52 | yes | yes | 2.5 | 2.0 |
| 5a7b | 31 | 10 | 1.94 | 3.69 | 1.94 | 1.64 | -7.01 | -6.78 | -10.26 | -10.26 | yes | yes | 10.0 | 8.9 |
| 3k5v | 27 | 5 | 3.34 | 3.33 | 2.14 | 2.24 | -9.82 | -9.84 | -11.43 | -11.43 | yes | yes | 4.4 | 3.1 |
| 3o9i | 42 | 16 | 3.91 | 3.94 | 1.93 | 0.81 | -9.35 | -9.38 | -17.74 | -17.75 | no | no | 26.5 | 23.2 |
| 3dx1 | 9 | 4 | 2.07 | 2.05 | 2.07 | 2.03 | -5.84 | -5.85 | -5.81 | -5.82 | yes | yes | 1.3 | 0.7 |
| 2vvn | 16 | 6 | 2.96 | 2.93 | 1.14 | 1.06 | -7.30 | -7.48 | -8.50 | -8.49 | yes | yes | 2.0 | 1.6 |
| 4eo8 | 27 | 7 | 1.01 | 1.01 | 0.97 | 1.01 | -7.04 | -7.06 | -10.15 | -10.15 | yes | yes | 5.9 | 4.0 |
| 3f3e | 9 | 6 | 0.78 | 0.56 | 0.57 | 0.56 | -6.67 | -6.73 | -7.67 | -7.67 | yes | yes | 1.1 | 0.6 |
| 3udh | 14 | 0 | 0.37 | 0.38 | 0.37 | 0.38 | -7.51 | -7.52 | -7.11 | -7.12 | yes | yes | 0.9 | 0.5 |
| 1r5y | 13 | 2 | 0.17 | 0.17 | 0.17 | 0.17 | -7.75 | -7.76 | -7.34 | -7.35 | yes | yes | 1.2 | 0.7 |
| 4gkm | 20 | 5 | 4.29 | 4.33 | 0.76 | 0.65 | -8.65 | -8.66 | -8.99 | -8.98 | yes | yes | 2.4 | 1.6 |
| 2hb1 | 14 | 4 | 4.19 | 4.20 | 1.26 | 1.19 | -6.33 | -6.30 | -6.87 | -6.87 | yes | yes | 1.6 | 1.1 |
| 5tmn | 32 | 17 | 2.20 | 2.75 | 2.20 | 1.78 | -8.34 | -8.68 | -12.15 | -12.14 | no | yes | 11.9 | 10.2 |
| 3pxf | 21 | 3 | 0.87 | 0.81 | 0.87 | 0.77 | -10.63 | -10.62 | -11.82 | -11.82 | yes | yes | 2.0 | 1.5 |
| 4jia | 31 | 7 | 1.31 | 1.32 | 0.97 | 0.95 | -9.78 | -9.82 | -11.97 | -11.97 | yes | yes | 4.7 | 4.6 |
| 3acw | 21 | 3 | 1.27 | 1.32 | 1.27 | 1.32 | -10.88 | -10.90 | -10.96 | -10.95 | yes | yes | 2.4 | 1.6 |
| 3u9q | 12 | 9 | 1.16 | 1.34 | 1.16 | 1.34 | -5.96 | -6.07 | -7.56 | -7.56 | yes | yes | 2.2 | 1.1 |
| 1owh | 24 | 5 | 3.17 | 3.16 | 0.87 | 0.88 | -8.62 | -8.61 | -9.30 | -9.30 | yes | yes | 3.7 | 2.8 |
| 3b27 | 16 | 4 | 0.48 | 0.44 | 0.48 | 0.44 | -7.18 | -7.19 | -7.67 | -7.67 | yes | yes | 1.6 | 1.0 |
| 1p1n | 15 | 5 | 0.99 | 1.00 | 0.99 | 1.00 | -8.33 | -8.29 | -8.88 | -8.87 | yes | yes | 1.6 | 1.3 |
| 3fv1 | 21 | 8 | 0.62 | 0.64 | 0.62 | 0.64 | -7.32 | -7.44 | -13.03 | -13.03 | yes | yes | 3.1 | 2.5 |
| 2y5h | 31 | 9 | 2.44 | 2.47 | 1.01 | 0.95 | -8.45 | -8.30 | -10.31 | -10.32 | yes | yes | 6.2 | 5.7 |
| 1uto | 9 | 3 | 1.04 | 1.04 | 1.01 | 1.04 | -4.53 | -4.54 | -4.79 | -4.79 | yes | yes | 1.1 | 0.5 |
| 4j21 | 20 | 3 | 0.24 | 0.26 | 0.24 | 0.26 | -11.84 | -11.85 | -12.69 | -12.69 | yes | yes | 1.9 | 1.4 |
| 3ozs | 33 | 10 | 1.03 | 1.00 | 1.03 | 1.00 | -9.82 | -9.76 | -14.61 | -14.62 | yes | yes | 11.3 | 12.1 |
| 3n76 | 22 | 8 | 0.66 | 0.64 | 0.66 | 0.64 | -11.00 | -11.00 | -14.14 | -14.14 | no | no | 3.7 | 3.0 |
| 3nq9 | 10 | 7 | 2.74 | 3.67 | 2.05 | 1.99 | -4.41 | -4.35 | -3.74 | -3.74 | yes | no | 1.4 | 0.7 |
| 1o0h | 27 | 9 | 2.97 | 2.84 | 1.25 | 1.37 | -6.48 | -6.36 | -7.68 | -7.68 | yes | yes | 6.7 | 5.9 |

## PoseBusters failures (top pose)

- `3o9i` vina: internal_steric_clash,internal_energy
- `3o9i` rustdock: internal_steric_clash,internal_energy
- `5tmn` vina: internal_steric_clash,internal_energy
- `3n76` vina: double_bond_flatness
- `3n76` rustdock: double_bond_flatness
- `3nq9` rustdock: internal_energy
