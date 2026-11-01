# Docking parity: rustdock-vina vs AutoDock Vina 1.2.7 (vina scoring)

CASF-2016 core set, subset `full` (234 complexes). Identical inputs for both: Meeko-prepared ligand from an RDKit conformer, Open Babel receptor, box = crystal ligand extent + 10 A, exhaustiveness 8, 9 modes, seed 42, 16 CPUs. Success = top pose within 2 A heavy-atom RMSD (symmetry-aware) of the crystal ligand.

Flexible side chains: up to 3 residues whose side chains have a heavy atom within 4 A of the crystal ligand, closest first, prepared with Meeko's `mk_prepare_receptor.py --flexres` (rigid part and flexible residues) from the same protein. The box covers the crystal ligand and the side chains plus 10 A. Complexes Meeko cannot prepare are left out for both programs (`flex_excluded.tsv` in the work directory). PoseBusters checks the top ligand pose against the protein with its crystal side chains.

| | AutoDock Vina | rustdock-vina |
|---|---|---|
| Docked without error | 234/234 | 234/234 |
| Top-1 success (< 2 A) | 101/234 | 104/234 |
| Any-mode success (< 2 A) | 172/234 | 166/234 |
| Median top-1 RMSD (A) | 2.72 | 2.63 |
| Top-1 PoseBusters-valid | 163/234 | 162/234 |
| Median wall time (s) | 15.5 | 10.5 |

Crystal-pose rescoring, rustdock minus Vina:

- intermolecular energy (scoring function alone): unavailable
- affinity (with torsional normalization): unavailable

Output round trip: Vina read rustdock's top pose from its output file in no complexes (not run).

| Code | Heavy | Rot | Vina RMSD | Rust RMSD | Vina best | Rust best | Vina aff | Rust aff | Xtal inter Vina | Xtal inter Rust | Vina PB | Rust PB | Vina s | Rust s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1a30 | 26 | 15 | 7.31 | 7.68 | 4.41 | 3.91 | -7.17 | -7.24 | - | - | yes | yes | 22.8 | 16.1 |
| 1bzc | 29 | 8 | 2.96 | 2.82 | 1.09 | 1.70 | -9.59 | -9.56 | - | - | no | no | 23.6 | 15.4 |
| 1c5z | 9 | 1 | 0.25 | 0.29 | 0.25 | 0.29 | -5.27 | -5.33 | - | - | no | yes | 3.6 | 1.8 |
| 1eby | 48 | 17 | 4.30 | 4.31 | 0.98 | 0.95 | -11.24 | -11.14 | - | - | yes | yes | 105.2 | 63.7 |
| 1g2k | 45 | 14 | 9.16 | 9.14 | 2.23 | 1.33 | -8.53 | -8.48 | - | - | yes | yes | 143.9 | 58.7 |
| 1k1i | 36 | 8 | 6.47 | 6.48 | 2.79 | 2.73 | -8.69 | -8.69 | - | - | no | no | 23.2 | 18.0 |
| 1lpg | 41 | 13 | 7.85 | 7.87 | 7.85 | 4.47 | -9.25 | -9.40 | - | - | yes | yes | 43.5 | 32.9 |
| 1mq6 | 36 | 10 | 2.81 | 8.32 | 1.30 | 1.33 | -8.78 | -8.84 | - | - | no | yes | 27.9 | 21.7 |
| 1nc1 | 20 | 7 | 0.83 | 0.84 | 0.83 | 0.84 | -8.00 | -8.02 | - | - | yes | yes | 13.3 | 7.7 |
| 1nc3 | 19 | 6 | 0.57 | 0.61 | 0.57 | 0.61 | -8.98 | -9.06 | - | - | yes | yes | 10.0 | 6.5 |
| 1nvq | 36 | 6 | 0.32 | 0.31 | 0.31 | 0.28 | -13.60 | -13.67 | - | - | yes | yes | 15.5 | 17.7 |
| 1o0h | 27 | 9 | 11.39 | 7.99 | 5.52 | 5.44 | -6.86 | -6.26 | - | - | no | no | 29.1 | 21.3 |
| 1o3f | 25 | 4 | 1.49 | 1.46 | 1.49 | 0.70 | -8.58 | -8.58 | - | - | yes | yes | 9.9 | 7.1 |
| 1owh | 24 | 5 | 2.50 | 2.37 | 1.01 | 1.00 | -8.43 | -8.34 | - | - | no | no | 13.3 | 8.7 |
| 1p1n | 15 | 5 | 0.99 | 1.00 | 0.99 | 1.00 | -8.04 | -8.03 | - | - | yes | yes | 12.6 | 6.9 |
| 1p1q | 13 | 5 | 4.67 | 4.67 | 0.81 | 0.72 | -7.87 | -7.86 | - | - | yes | yes | 11.9 | 7.2 |
| 1ps3 | 16 | 5 | 0.67 | 0.98 | 0.67 | 0.98 | -7.20 | -7.42 | - | - | yes | no | 9.8 | 6.6 |
| 1pxn | 22 | 7 | 1.29 | 1.30 | 1.08 | 1.29 | -8.87 | -8.56 | - | - | yes | yes | 13.2 | 8.4 |
| 1q8t | 18 | 5 | 7.65 | 7.64 | 1.74 | 1.17 | -8.05 | -8.10 | - | - | yes | yes | 8.0 | 5.4 |
| 1q8u | 22 | 4 | 1.85 | 1.88 | 1.85 | 1.88 | -6.56 | -6.56 | - | - | yes | yes | 22.9 | 9.2 |
| 1qf1 | 26 | 14 | 2.60 | 2.62 | 0.89 | 0.98 | -8.00 | -8.10 | - | - | yes | yes | 37.4 | 26.6 |
| 1qkt | 20 | 3 | 1.25 | 1.21 | 1.25 | 0.61 | -10.47 | -10.49 | - | - | no | no | 12.1 | 8.4 |
| 1r5y | 13 | 2 | 0.54 | 0.46 | 0.54 | 0.46 | -7.48 | -7.49 | - | - | yes | yes | 6.5 | 3.8 |
| 1s38 | 13 | 2 | 4.95 | 3.50 | 0.55 | 0.54 | -7.83 | -7.68 | - | - | no | no | 5.7 | 3.4 |
| 1sqa | 31 | 7 | 2.63 | 2.64 | 1.06 | 1.09 | -9.28 | -9.27 | - | - | yes | yes | 19.1 | 14.1 |
| 1syi | 17 | 4 | 0.23 | 0.24 | 0.23 | 0.24 | -8.84 | -8.83 | - | - | yes | yes | 13.7 | 10.1 |
| 1u1b | 51 | 17 | 11.06 | 11.02 | 3.86 | 7.81 | -8.08 | -8.09 | - | - | yes | no | 78.6 | 55.2 |
| 1uto | 9 | 3 | 2.07 | 1.65 | 1.08 | 1.10 | -4.67 | -4.61 | - | - | no | no | 3.8 | 1.9 |
| 1vso | 21 | 7 | 4.96 | 1.41 | 1.70 | 1.41 | -6.62 | -7.05 | - | - | yes | no | 16.2 | 9.9 |
| 1w4o | 21 | 6 | 1.83 | 1.83 | 1.83 | 1.83 | -6.04 | -6.05 | - | - | yes | yes | 15.1 | 11.1 |
| 1y6r | 20 | 7 | 0.99 | 0.99 | 0.99 | 0.99 | -9.32 | -9.32 | - | - | no | no | 13.8 | 9.5 |
| 1yc1 | 27 | 7 | 0.83 | 0.85 | 0.83 | 0.85 | -9.93 | -9.98 | - | - | yes | yes | 17.0 | 11.6 |
| 1ydr | 20 | 3 | 3.45 | 3.42 | 1.81 | 2.17 | -7.92 | -7.95 | - | - | yes | yes | 9.8 | 9.3 |
| 1ydt | 27 | 8 | 4.59 | 9.64 | 3.02 | 4.12 | -8.54 | -8.63 | - | - | no | no | 15.6 | 12.0 |
| 1z6e | 38 | 9 | 1.39 | 1.42 | 1.39 | 1.42 | -10.92 | -10.92 | - | - | yes | yes | 26.8 | 21.6 |
| 1z95 | 29 | 7 | 0.67 | 0.67 | 0.41 | 0.67 | -10.38 | -10.43 | - | - | yes | yes | 30.5 | 22.4 |
| 1z9g | 17 | 7 | 3.68 | 3.63 | 3.65 | 0.84 | -6.23 | -6.22 | - | - | yes | yes | 16.9 | 10.7 |
| 2al5 | 15 | 4 | 0.51 | 0.46 | 0.51 | 0.46 | -8.42 | -8.39 | - | - | yes | yes | 13.5 | 8.2 |
| 2br1 | 29 | 10 | 7.74 | 7.42 | 1.36 | 1.54 | -8.16 | -8.03 | - | - | yes | yes | 16.5 | 12.0 |
| 2c3i | 23 | 6 | 5.54 | 5.54 | 2.55 | 1.82 | -9.50 | -9.45 | - | - | yes | yes | 15.9 | 9.8 |
| 2cbv | 12 | 4 | 3.16 | 3.16 | 1.49 | 1.47 | -7.25 | -7.25 | - | - | yes | yes | 9.3 | 6.3 |
| 2cet | 22 | 8 | 3.28 | 8.17 | 3.28 | 3.30 | -8.95 | -8.54 | - | - | no | no | 17.4 | 11.2 |
| 2fvd | 30 | 9 | 9.68 | 9.79 | 3.06 | 3.17 | -8.50 | -8.43 | - | - | no | no | 29.3 | 21.0 |
| 2fxs | 28 | 12 | 5.89 | 7.56 | 3.44 | 4.16 | -7.35 | -7.27 | - | - | no | no | 31.7 | 19.4 |
| 2hb1 | 14 | 4 | 3.85 | 3.83 | 0.97 | 1.26 | -6.33 | -6.32 | - | - | yes | yes | 15.7 | 8.4 |
| 2iwx | 24 | 2 | 0.69 | 0.71 | 0.69 | 0.71 | -5.27 | -5.21 | - | - | no | yes | 46.3 | 17.9 |
| 2j78 | 13 | 6 | 9.18 | 9.21 | 3.26 | 3.24 | -6.25 | -6.23 | - | - | yes | yes | 12.2 | 6.6 |
| 2p15 | 32 | 5 | 0.43 | 0.60 | 0.43 | 0.55 | -13.60 | -13.65 | - | - | no | no | 25.9 | 18.9 |
| 2p4y | 40 | 11 | 0.55 | 0.61 | 0.55 | 0.61 | -10.76 | -10.77 | - | - | yes | yes | 31.3 | 25.0 |
| 2pog | 21 | 3 | 3.42 | 3.44 | 0.69 | 0.71 | -10.68 | -10.63 | - | - | no | no | 14.4 | 9.7 |
| 2qbp | 37 | 10 | 1.10 | 1.73 | 1.10 | 1.73 | -8.77 | -8.34 | - | - | yes | yes | 45.0 | 30.6 |
| 2qbq | 31 | 11 | 1.01 | 0.89 | 1.01 | 0.89 | -8.87 | -8.81 | - | - | no | no | 28.2 | 21.6 |
| 2qbr | 28 | 8 | 0.96 | 0.96 | 0.96 | 0.96 | -8.59 | -8.79 | - | - | no | no | 26.6 | 17.2 |
| 2qe4 | 24 | 6 | 0.51 | 0.49 | 0.50 | 0.41 | -9.55 | -9.56 | - | - | yes | yes | 18.6 | 12.3 |
| 2qnq | 41 | 10 | 0.56 | 0.62 | 0.56 | 0.62 | -10.15 | -10.25 | - | - | yes | yes | 35.4 | 26.3 |
| 2v00 | 16 | 4 | 1.00 | 1.01 | 1.00 | 1.01 | -7.55 | -7.54 | - | - | yes | yes | 7.3 | 4.5 |
| 2v7a | 35 | 8 | 2.23 | 2.17 | 2.05 | 2.17 | -9.28 | -9.26 | - | - | no | no | 25.3 | 17.6 |
| 2vkm | 47 | 20 | 2.37 | 8.09 | 2.37 | 7.31 | -11.63 | -9.69 | - | - | yes | yes | 99.0 | 69.1 |
| 2vw5 | 36 | 10 | 5.61 | 5.73 | 2.87 | 3.09 | -6.28 | -5.97 | - | - | no | no | 118.6 | 56.8 |
| 2w4x | 18 | 8 | 1.86 | 1.87 | 1.86 | 1.87 | -6.22 | -6.25 | - | - | yes | yes | 14.6 | 10.1 |
| 2wbg | 22 | 12 | 8.98 | 8.41 | 3.28 | 3.23 | -7.69 | -7.52 | - | - | yes | yes | 22.0 | 13.7 |
| 2wca | 27 | 10 | 5.07 | 5.07 | 5.07 | 5.07 | -7.86 | -7.80 | - | - | no | no | 21.4 | 15.1 |
| 2weg | 11 | 2 | 2.34 | 2.00 | 2.28 | 1.19 | -6.16 | -6.12 | - | - | yes | yes | 5.3 | 3.3 |
| 2wer | 25 | 3 | 4.72 | 6.55 | 1.77 | 2.59 | -5.38 | -5.12 | - | - | no | no | 50.2 | 19.2 |
| 2wtv | 34 | 4 | 5.97 | 6.02 | 3.82 | 3.96 | -9.45 | -9.41 | - | - | no | no | 18.7 | 14.9 |
| 2wvt | 13 | 7 | 1.51 | 1.48 | 1.46 | 1.48 | -6.65 | -6.66 | - | - | yes | yes | 10.6 | 6.6 |
| 2x00 | 37 | 7 | 4.17 | 4.59 | 1.51 | 4.02 | -7.97 | -7.70 | - | - | no | no | 114.6 | 46.2 |
| 2xb8 | 22 | 8 | 1.54 | 1.55 | 1.54 | 1.39 | -10.05 | -10.05 | - | - | yes | yes | 13.6 | 9.6 |
| 2xdl | 16 | 8 | 0.96 | 0.94 | 0.96 | 0.94 | -6.24 | -6.24 | - | - | yes | yes | 8.4 | 5.3 |
| 2xii | 28 | 7 | 6.86 | 6.29 | 3.64 | 3.67 | -8.61 | -8.60 | - | - | yes | yes | 17.6 | 13.4 |
| 2xj7 | 16 | 5 | 3.01 | 3.02 | 3.01 | 3.02 | -6.78 | -6.81 | - | - | yes | yes | 7.6 | 5.4 |
| 2xnb | 27 | 6 | 0.92 | 0.84 | 0.92 | 0.84 | -10.36 | -10.43 | - | - | yes | yes | 13.6 | 9.2 |
| 2xys | 25 | 0 | 4.76 | 4.76 | 0.59 | 0.58 | -9.35 | -9.29 | - | - | yes | yes | 43.8 | 19.5 |
| 2yfe | 25 | 11 | 5.40 | 5.67 | 1.16 | 1.14 | -8.19 | -8.29 | - | - | yes | no | 12.7 | 9.8 |
| 2yge | 40 | 12 | 7.83 | 3.17 | 2.23 | 2.06 | -5.44 | -5.61 | - | - | no | no | 120.6 | 53.1 |
| 2yki | 34 | 3 | 2.51 | 2.51 | 2.51 | 2.51 | -14.74 | -14.74 | - | - | yes | yes | 21.1 | 15.9 |
| 2ymd | 13 | 4 | 0.29 | 0.27 | 0.29 | 0.27 | -7.25 | -7.25 | - | - | yes | yes | 17.4 | 14.1 |
| 2zb1 | 26 | 7 | 1.38 | 1.38 | 1.38 | 1.38 | -10.99 | -10.98 | - | - | yes | yes | 10.7 | 7.6 |
| 2zcq | 25 | 9 | 5.50 | 5.44 | 3.19 | 2.10 | -9.50 | -9.71 | - | - | yes | yes | 19.4 | 11.2 |
| 2zcr | 28 | 12 | 3.96 | 3.97 | 1.45 | 1.40 | -10.17 | -10.24 | - | - | yes | yes | 20.0 | 12.4 |
| 2zy1 | 24 | 9 | 6.30 | 3.24 | 2.94 | 2.12 | -9.47 | -9.47 | - | - | yes | yes | 25.7 | 16.1 |
| 3acw | 21 | 3 | 1.28 | 1.28 | 1.28 | 1.28 | -10.77 | -10.86 | - | - | yes | yes | 10.9 | 7.5 |
| 3ag9 | 67 | 36 | 9.31 | 9.55 | 8.39 | 7.44 | -7.67 | -7.82 | - | - | yes | yes | 353.1 | 219.8 |
| 3arp | 34 | 15 | 7.26 | 7.21 | 6.74 | 5.65 | -9.53 | -9.35 | - | - | no | no | 58.0 | 38.8 |
| 3arq | 36 | 8 | 9.55 | 9.51 | 2.23 | 2.19 | -11.33 | -11.34 | - | - | no | no | 34.8 | 31.5 |
| 3aru | 20 | 8 | 11.88 | 11.84 | 2.12 | 1.61 | -6.63 | -6.62 | - | - | yes | yes | 15.8 | 11.3 |
| 3arv | 25 | 1 | 9.36 | 9.34 | 3.89 | 4.60 | -11.41 | -11.38 | - | - | yes | yes | 20.0 | 18.9 |
| 3ary | 17 | 2 | 5.68 | 5.62 | 5.28 | 5.62 | -7.25 | -7.23 | - | - | yes | yes | 12.8 | 9.9 |
| 3b1m | 38 | 12 | 2.22 | 2.22 | 0.65 | 2.22 | -12.10 | -12.09 | - | - | no | no | 24.7 | 21.6 |
| 3b27 | 16 | 4 | 0.58 | 0.46 | 0.58 | 0.46 | -7.04 | -7.04 | - | - | yes | yes | 7.2 | 4.8 |
| 3b5r | 28 | 7 | 0.57 | 0.58 | 0.55 | 0.58 | -10.74 | -10.80 | - | - | yes | yes | 28.4 | 20.6 |
| 3b65 | 25 | 7 | 0.53 | 0.55 | 0.29 | 0.36 | -10.09 | -10.03 | - | - | yes | yes | 23.6 | 19.9 |
| 3b68 | 31 | 10 | 0.72 | 0.72 | 0.72 | 0.72 | -11.21 | -11.31 | - | - | yes | yes | 30.5 | 14.4 |
| 3bgz | 24 | 3 | 0.79 | 0.63 | 0.79 | 0.63 | -11.38 | -11.32 | - | - | no | no | 14.7 | 10.0 |
| 3coy | 31 | 10 | 4.05 | 4.05 | 1.39 | 0.42 | -9.61 | -9.60 | - | - | no | yes | 31.7 | 23.1 |
| 3coz | 30 | 13 | 4.84 | 4.86 | 1.49 | 0.93 | -9.63 | -9.66 | - | - | no | no | 29.8 | 21.9 |
| 3d4z | 14 | 5 | 1.00 | 0.99 | 1.00 | 0.99 | -7.09 | -7.07 | - | - | yes | yes | 7.9 | 4.8 |
| 3d6q | 22 | 5 | 6.76 | 5.07 | 4.02 | 3.92 | -6.30 | -6.28 | - | - | no | no | 15.6 | 10.6 |
| 3dd0 | 16 | 5 | 1.08 | 1.07 | 1.08 | 1.07 | -6.57 | -6.58 | - | - | yes | yes | 7.8 | 4.9 |
| 3dx2 | 12 | 6 | 3.53 | 3.56 | 3.52 | 3.43 | -5.07 | -5.06 | - | - | yes | yes | 9.5 | 5.5 |
| 3dxg | 21 | 6 | 2.81 | 2.80 | 2.81 | 2.80 | -6.55 | -6.55 | - | - | yes | yes | 13.5 | 10.4 |
| 3e92 | 26 | 7 | 0.24 | 0.23 | 0.24 | 0.23 | -12.74 | -12.78 | - | - | yes | yes | 13.8 | 10.0 |
| 3e93 | 37 | 6 | 3.30 | 3.32 | 0.71 | 0.70 | -13.48 | -13.49 | - | - | yes | yes | 34.4 | 22.4 |
| 3ebp | 28 | 6 | 1.76 | 1.71 | 1.76 | 1.71 | -9.54 | -9.46 | - | - | no | no | 25.8 | 21.5 |
| 3ehy | 17 | 7 | 1.82 | 1.79 | 1.82 | 1.49 | -7.17 | -7.15 | - | - | yes | yes | 11.0 | 7.5 |
| 3fcq | 14 | 4 | 9.21 | 9.18 | 1.58 | 0.69 | -5.60 | -5.59 | - | - | yes | yes | 8.1 | 4.8 |
| 3fur | 31 | 5 | 8.37 | 0.78 | 1.62 | 0.78 | -10.98 | -11.00 | - | - | no | yes | 31.8 | 27.3 |
| 3fv1 | 21 | 8 | 5.08 | 5.07 | 0.71 | 0.70 | -7.46 | -7.50 | - | - | no | no | 15.0 | 9.0 |
| 3fv2 | 20 | 7 | 0.93 | 0.94 | 0.93 | 0.94 | -8.33 | -8.40 | - | - | no | no | 18.3 | 11.1 |
| 3g0w | 24 | 3 | 0.73 | 0.80 | 0.73 | 0.80 | -9.69 | -9.69 | - | - | yes | yes | 19.1 | 12.6 |
| 3g2n | 20 | 7 | 1.67 | 1.66 | 1.67 | 1.66 | -9.02 | -8.99 | - | - | yes | yes | 14.0 | 9.6 |
| 3g2z | 13 | 2 | 5.30 | 5.34 | 3.07 | 3.77 | -6.53 | -6.55 | - | - | yes | yes | 4.8 | 2.6 |
| 3g31 | 17 | 3 | 4.86 | 3.60 | 1.20 | 3.35 | -4.64 | -4.41 | - | - | no | no | 11.5 | 7.7 |
| 3gbb | 18 | 5 | 1.19 | 1.18 | 1.19 | 1.18 | -8.06 | -8.02 | - | - | yes | yes | 9.2 | 5.9 |
| 3gc5 | 20 | 6 | 5.46 | 4.00 | 1.03 | 0.56 | -8.57 | -8.73 | - | - | no | no | 9.5 | 6.6 |
| 3ge7 | 26 | 8 | 0.99 | 0.98 | 0.99 | 0.98 | -10.11 | -10.15 | - | - | yes | yes | 14.8 | 11.0 |
| 3gnw | 44 | 9 | 1.65 | 7.52 | 1.65 | 5.23 | -8.70 | -7.52 | - | - | no | no | 70.2 | 29.3 |
| 3gr2 | 14 | 4 | 7.06 | 6.92 | 3.64 | 1.72 | -6.07 | -6.03 | - | - | yes | yes | 7.4 | 5.0 |
| 3gv9 | 12 | 3 | 4.56 | 4.68 | 1.62 | 1.19 | -5.28 | -5.38 | - | - | yes | yes | 5.7 | 3.3 |
| 3gy4 | 10 | 2 | 0.58 | 0.60 | 0.58 | 0.60 | -5.85 | -5.84 | - | - | yes | yes | 4.5 | 2.2 |
| 3ivg | 29 | 8 | 1.56 | 1.54 | 1.39 | 1.38 | -9.93 | -9.91 | - | - | yes | yes | 22.9 | 13.9 |
| 3jvr | 23 | 4 | 16.80 | 17.08 | 0.78 | 0.77 | -7.25 | -7.38 | - | - | yes | yes | 15.5 | 13.1 |
| 3jvs | 30 | 3 | 9.89 | 9.90 | 2.40 | 2.45 | -8.08 | -8.03 | - | - | yes | yes | 24.7 | 19.0 |
| 3jya | 16 | 0 | 4.62 | 4.60 | 2.97 | 0.43 | -7.97 | -8.01 | - | - | yes | yes | 4.9 | 3.4 |
| 3kgp | 11 | 3 | 1.15 | 1.16 | 1.06 | 1.00 | -5.99 | -5.98 | - | - | yes | yes | 5.2 | 2.7 |
| 3kr8 | 21 | 1 | 0.41 | 0.44 | 0.41 | 0.44 | -11.49 | -11.50 | - | - | yes | yes | 9.8 | 7.2 |
| 3l7b | 19 | 6 | 1.17 | 1.13 | 1.17 | 1.13 | -7.57 | -7.56 | - | - | yes | yes | 12.3 | 8.6 |
| 3lka | 12 | 4 | 5.52 | 5.53 | 1.93 | 1.20 | -6.60 | -6.62 | - | - | yes | yes | 7.9 | 4.6 |
| 3mss | 21 | 8 | 2.46 | 8.23 | 2.08 | 1.50 | -6.92 | -6.86 | - | - | yes | yes | 9.9 | 7.0 |
| 3myg | 30 | 13 | 12.55 | 12.55 | 1.40 | 10.65 | -8.25 | -8.12 | - | - | no | no | 28.2 | 19.4 |
| 3n76 | 22 | 8 | 0.72 | 0.68 | 0.72 | 0.68 | -10.65 | -10.68 | - | - | no | no | 21.0 | 15.0 |
| 3n7a | 12 | 4 | 0.99 | 1.00 | 0.99 | 1.00 | -6.64 | -6.65 | - | - | yes | yes | 6.1 | 3.8 |
| 3n86 | 23 | 7 | 0.61 | 0.62 | 0.61 | 0.62 | -11.36 | -11.38 | - | - | yes | yes | 12.8 | 9.7 |
| 3nq9 | 10 | 7 | 6.80 | 7.12 | 2.50 | 2.06 | -5.49 | -5.57 | - | - | yes | yes | 8.7 | 4.0 |
| 3nw9 | 38 | 11 | 0.49 | 0.53 | 0.49 | 0.53 | -11.72 | -11.74 | - | - | yes | yes | 36.1 | 26.4 |
| 3nx7 | 20 | 10 | 1.14 | 1.21 | 1.14 | 1.00 | -6.85 | -6.80 | - | - | yes | yes | 16.6 | 11.7 |
| 3o9i | 42 | 16 | 5.16 | 6.85 | 0.77 | 1.42 | -9.28 | -8.98 | - | - | yes | yes | 57.9 | 61.7 |
| 3oe4 | 33 | 10 | 0.65 | 0.65 | 0.65 | 0.65 | -9.97 | -9.91 | - | - | no | no | 33.6 | 28.1 |
| 3oe5 | 31 | 11 | 1.74 | 1.56 | 1.74 | 1.56 | -8.75 | -8.78 | - | - | yes | no | 34.2 | 23.1 |
| 3ozs | 33 | 10 | 1.11 | 11.46 | 1.11 | 1.12 | -9.44 | -9.06 | - | - | yes | no | 33.3 | 24.6 |
| 3ozt | 31 | 10 | 0.42 | 0.58 | 0.42 | 0.58 | -10.10 | -9.77 | - | - | no | yes | 28.4 | 18.4 |
| 3p5o | 30 | 8 | 5.12 | 5.12 | 5.08 | 5.08 | -9.27 | -9.26 | - | - | no | no | 17.6 | 13.7 |
| 3prs | 50 | 23 | 4.48 | 4.50 | 4.48 | 4.50 | -9.63 | -9.66 | - | - | yes | yes | 103.7 | 69.8 |
| 3pww | 49 | 13 | 2.45 | 2.46 | 1.25 | 1.61 | -10.86 | -10.88 | - | - | yes | yes | 82.1 | 47.7 |
| 3pxf | 21 | 3 | 1.07 | 1.08 | 1.07 | 1.08 | -10.69 | -10.64 | - | - | yes | yes | 9.2 | 6.1 |
| 3qgy | 30 | 6 | 8.94 | 8.94 | 2.21 | 2.19 | -9.67 | -9.67 | - | - | yes | yes | 16.8 | 12.0 |
| 3qqs | 19 | 4 | 4.47 | 4.90 | 4.41 | 4.37 | -8.92 | -7.95 | - | - | yes | no | 21.4 | 13.2 |
| 3r88 | 14 | 6 | 14.42 | 13.92 | 2.15 | 5.49 | -6.24 | -6.13 | - | - | yes | yes | 15.9 | 8.8 |
| 3rlr | 23 | 5 | 0.44 | 0.46 | 0.44 | 0.46 | -8.71 | -8.68 | - | - | yes | yes | 8.3 | 6.0 |
| 3rr4 | 18 | 4 | 4.79 | 4.80 | 0.73 | 3.31 | -8.65 | -8.69 | - | - | no | no | 9.6 | 6.9 |
| 3rsx | 16 | 2 | 5.43 | 5.41 | 0.52 | 3.55 | -8.52 | -8.53 | - | - | yes | yes | 9.2 | 6.2 |
| 3ryj | 18 | 4 | 3.44 | 3.42 | 0.97 | 0.96 | -7.38 | -7.39 | - | - | yes | yes | 8.9 | 5.8 |
| 3syr | 20 | 6 | 0.63 | 0.63 | 0.63 | 0.63 | -9.84 | -9.88 | - | - | yes | yes | 11.6 | 7.8 |
| 3tsk | 40 | 15 | 5.23 | 5.31 | 1.66 | 3.79 | -9.79 | -9.70 | - | - | yes | yes | 67.2 | 58.0 |
| 3twp | 10 | 2 | 6.65 | 6.74 | 2.12 | 2.00 | -5.75 | -5.79 | - | - | yes | yes | 12.3 | 6.8 |
| 3u5j | 22 | 2 | 5.07 | 5.06 | 4.64 | 4.70 | -8.80 | -8.76 | - | - | yes | yes | 9.4 | 6.8 |
| 3u8n | 17 | 4 | 9.29 | 0.52 | 1.65 | 0.52 | -4.54 | -5.32 | - | - | yes | yes | 48.5 | 19.6 |
| 3u9q | 12 | 9 | 1.26 | 1.16 | 1.26 | 1.16 | -6.09 | -6.04 | - | - | yes | yes | 12.7 | 8.1 |
| 3udh | 14 | 0 | 0.42 | 0.43 | 0.42 | 0.43 | -7.56 | -7.55 | - | - | yes | yes | 5.8 | 4.0 |
| 3ueu | 14 | 11 | 4.46 | 4.47 | 1.80 | 1.92 | -6.38 | -6.33 | - | - | no | no | 11.0 | 5.9 |
| 3uev | 16 | 13 | 5.05 | 5.17 | 1.38 | 2.24 | -6.44 | -6.70 | - | - | no | no | 18.0 | 9.5 |
| 3uew | 18 | 15 | 11.34 | 2.43 | 2.19 | 2.43 | -6.35 | -6.35 | - | - | yes | yes | 21.5 | 11.0 |
| 3uex | 20 | 17 | 4.33 | 4.96 | 2.69 | 2.53 | -7.20 | -7.29 | - | - | yes | no | 20.7 | 10.3 |
| 3ui7 | 24 | 6 | 3.46 | 3.63 | 3.46 | 0.93 | -9.01 | -9.03 | - | - | no | no | 13.9 | 11.2 |
| 3uo4 | 29 | 6 | 5.87 | 5.86 | 1.96 | 1.74 | -9.62 | -9.68 | - | - | yes | yes | 30.3 | 18.4 |
| 3up2 | 28 | 6 | 2.59 | 2.61 | 0.94 | 0.85 | -9.61 | -9.73 | - | - | yes | yes | 28.5 | 19.3 |
| 3uuo | 23 | 5 | 3.27 | 3.24 | 0.76 | 1.65 | -9.01 | -8.96 | - | - | no | no | 9.3 | 7.0 |
| 3wz8 | 33 | 9 | 5.92 | 5.87 | 4.16 | 5.87 | -8.95 | -8.90 | - | - | yes | yes | 21.5 | 19.2 |
| 3zdg | 13 | 9 | 0.57 | 0.56 | 0.57 | 0.56 | -5.55 | -5.54 | - | - | yes | yes | 16.9 | 12.5 |
| 3zsx | 34 | 10 | 4.63 | 4.67 | 2.17 | 2.15 | -7.43 | -7.46 | - | - | yes | yes | 36.6 | 24.1 |
| 3zt2 | 23 | 2 | 3.22 | 3.18 | 1.37 | 1.33 | -9.05 | -9.03 | - | - | no | no | 13.3 | 10.4 |
| 4agn | 24 | 9 | 1.83 | 1.84 | 1.83 | 1.84 | -5.62 | -5.59 | - | - | yes | yes | 13.1 | 9.0 |
| 4agp | 30 | 10 | 3.87 | 5.98 | 1.79 | 2.63 | -5.58 | -5.22 | - | - | no | no | 18.2 | 14.5 |
| 4agq | 30 | 10 | 1.95 | 1.93 | 1.95 | 1.63 | -6.12 | -6.13 | - | - | yes | yes | 19.1 | 14.8 |
| 4cig | 37 | 13 | 6.27 | 6.29 | 1.81 | 4.67 | -8.33 | -8.33 | - | - | no | no | 33.1 | 26.4 |
| 4cra | 38 | 9 | 10.14 | 10.21 | 1.58 | 1.61 | -9.63 | -10.23 | - | - | no | no | 30.2 | 25.4 |
| 4ddh | 15 | 4 | 2.02 | 2.08 | 1.61 | 1.76 | -6.68 | -6.56 | - | - | no | no | 9.0 | 5.2 |
| 4ddk | 12 | 1 | 4.20 | 4.77 | 1.00 | 2.09 | -6.07 | -6.03 | - | - | yes | yes | 4.7 | 3.1 |
| 4de2 | 24 | 7 | 0.55 | 1.51 | 0.55 | 1.35 | -8.17 | -8.07 | - | - | yes | no | 9.3 | 6.1 |
| 4de3 | 21 | 3 | 1.39 | 1.32 | 1.39 | 1.32 | -7.98 | -8.07 | - | - | yes | yes | 6.8 | 4.2 |
| 4djv | 28 | 6 | 1.04 | 1.04 | 1.04 | 1.04 | -9.62 | -9.65 | - | - | yes | yes | 14.8 | 11.2 |
| 4dld | 21 | 8 | 0.69 | 1.53 | 0.69 | 1.53 | -8.19 | -8.54 | - | - | yes | no | 10.2 | 5.7 |
| 4dli | 21 | 4 | 5.10 | 5.58 | 0.66 | 1.97 | -10.00 | -9.83 | - | - | no | yes | 12.4 | 9.1 |
| 4e5w | 28 | 5 | 1.06 | 1.08 | 1.06 | 1.08 | -9.69 | -9.68 | - | - | no | no | 19.9 | 14.1 |
| 4e6q | 25 | 3 | 2.83 | 2.66 | 0.97 | 0.65 | -9.44 | -9.18 | - | - | yes | yes | 10.3 | 8.4 |
| 4ea2 | 24 | 12 | 1.94 | 1.82 | 1.61 | 1.82 | -7.99 | -7.73 | - | - | no | no | 25.7 | 18.5 |
| 4eky | 24 | 8 | 0.77 | 0.76 | 0.77 | 0.76 | -9.15 | -9.15 | - | - | yes | yes | 23.0 | 13.2 |
| 4eo8 | 27 | 7 | 1.16 | 1.12 | 1.16 | 1.12 | -7.24 | -7.21 | - | - | yes | yes | 24.9 | 18.4 |
| 4eor | 28 | 7 | 1.42 | 1.42 | 1.42 | 1.42 | -9.14 | -9.17 | - | - | yes | yes | 17.5 | 12.2 |
| 4f09 | 19 | 2 | 3.41 | 3.41 | 0.41 | 0.35 | -7.76 | -7.76 | - | - | yes | yes | 5.9 | 4.0 |
| 4f2w | 20 | 7 | 0.35 | 0.35 | 0.35 | 0.35 | -8.46 | -8.41 | - | - | yes | yes | 10.5 | 7.6 |
| 4f3c | 23 | 10 | 0.81 | 0.81 | 0.81 | 0.81 | -9.56 | -9.39 | - | - | yes | yes | 20.1 | 14.1 |
| 4f9w | 25 | 5 | 0.50 | 0.47 | 0.50 | 0.47 | -11.78 | -11.86 | - | - | no | no | 13.2 | 9.3 |
| 4gfm | 18 | 2 | 6.31 | 6.32 | 0.47 | 1.93 | -7.50 | -7.49 | - | - | yes | yes | 9.2 | 5.8 |
| 4gid | 47 | 23 | 9.07 | 0.66 | 6.49 | 0.66 | -9.49 | -11.99 | - | - | no | yes | 91.7 | 74.5 |
| 4gkm | 20 | 5 | 5.76 | 5.75 | 5.33 | 4.86 | -8.50 | -8.51 | - | - | no | no | 26.4 | 17.7 |
| 4gr0 | 44 | 14 | 1.71 | 1.64 | 1.09 | 1.64 | -11.29 | -11.05 | - | - | yes | yes | 74.9 | 52.5 |
| 4hge | 25 | 4 | 5.88 | 1.63 | 0.47 | 1.63 | -8.53 | -9.17 | - | - | no | yes | 8.2 | 6.2 |
| 4ih5 | 17 | 3 | 1.73 | 1.76 | 1.73 | 1.76 | -6.53 | -6.55 | - | - | yes | yes | 9.9 | 6.5 |
| 4ih7 | 17 | 1 | 0.80 | 0.79 | 0.80 | 0.79 | -7.54 | -7.53 | - | - | yes | yes | 8.0 | 5.6 |
| 4ivb | 23 | 4 | 1.15 | 1.11 | 1.15 | 1.10 | -9.15 | -9.17 | - | - | yes | yes | 10.9 | 8.1 |
| 4ivc | 24 | 5 | 0.93 | 0.91 | 0.93 | 0.91 | -8.89 | -8.91 | - | - | yes | yes | 11.0 | 7.6 |
| 4ivd | 25 | 6 | 2.56 | 2.55 | 1.97 | 1.46 | -8.60 | -8.62 | - | - | yes | yes | 13.4 | 9.0 |
| 4j21 | 20 | 3 | 0.53 | 0.41 | 0.53 | 0.41 | -11.75 | -11.64 | - | - | yes | yes | 10.4 | 7.5 |
| 4j28 | 17 | 4 | 2.44 | 2.35 | 1.10 | 2.35 | -7.07 | -7.08 | - | - | yes | yes | 11.0 | 7.8 |
| 4j3l | 26 | 7 | 1.48 | 1.67 | 1.48 | 1.67 | -11.00 | -11.01 | - | - | no | yes | 18.9 | 14.6 |
| 4jfs | 17 | 7 | 1.12 | 1.12 | 1.12 | 1.12 | -6.76 | -6.82 | - | - | yes | yes | 9.7 | 6.2 |
| 4jia | 31 | 7 | 8.54 | 8.58 | 0.93 | 1.29 | -8.87 | -8.88 | - | - | yes | yes | 16.5 | 13.5 |
| 4jsz | 10 | 2 | 2.08 | 2.09 | 2.08 | 2.09 | -4.79 | -4.77 | - | - | no | no | 9.0 | 6.6 |
| 4jxs | 22 | 6 | 3.10 | 3.10 | 1.85 | 2.02 | -7.65 | -7.69 | - | - | yes | yes | 11.8 | 8.2 |
| 4k18 | 29 | 4 | 1.57 | 1.51 | 1.57 | 1.51 | -9.62 | -9.76 | - | - | no | no | 27.8 | 19.4 |
| 4k77 | 18 | 2 | 3.08 | 3.10 | 2.42 | 2.38 | -7.85 | -7.86 | - | - | yes | yes | 5.8 | 3.8 |
| 4kz6 | 15 | 5 | 4.86 | 4.85 | 1.15 | 1.20 | -5.43 | -5.39 | - | - | yes | yes | 8.2 | 5.1 |
| 4kzq | 18 | 2 | 0.82 | 0.90 | 0.82 | 0.90 | -10.52 | -10.53 | - | - | yes | yes | 11.5 | 8.5 |
| 4kzu | 18 | 1 | 0.29 | 0.31 | 0.29 | 0.31 | -10.80 | -10.80 | - | - | yes | yes | 11.1 | 7.9 |
| 4lzs | 15 | 7 | 0.44 | 0.40 | 0.44 | 0.40 | -6.26 | -6.23 | - | - | yes | yes | 7.0 | 4.3 |
| 4m0y | 22 | 3 | 7.51 | 7.53 | 5.05 | 5.98 | -8.83 | -8.80 | - | - | no | no | 14.0 | 9.4 |
| 4m0z | 24 | 5 | 0.44 | 0.42 | 0.44 | 0.42 | -10.52 | -10.57 | - | - | yes | yes | 15.3 | 9.8 |
| 4mgd | 19 | 4 | 0.62 | 0.63 | 0.62 | 0.63 | -8.83 | -8.85 | - | - | yes | yes | 14.5 | 9.1 |
| 4ogj | 37 | 11 | 7.66 | 7.64 | 4.53 | 4.61 | -7.93 | -7.79 | - | - | no | no | 36.0 | 26.5 |
| 4owm | 11 | 2 | 3.49 | 3.46 | 2.56 | 2.57 | -5.44 | -5.45 | - | - | yes | yes | 11.9 | 6.2 |
| 4pcs | 16 | 3 | 1.09 | 1.11 | 1.09 | 1.11 | -6.49 | -6.51 | - | - | yes | yes | 9.8 | 6.5 |
| 4qd6 | 35 | 8 | 4.92 | 4.92 | 1.85 | 1.33 | -8.14 | -8.10 | - | - | yes | yes | 28.3 | 21.8 |
| 4rfm | 36 | 6 | 4.67 | 4.10 | 0.88 | 0.77 | -11.08 | -11.24 | - | - | yes | yes | 15.4 | 15.9 |
| 4tmn | 36 | 16 | 9.15 | 9.12 | 1.40 | 1.32 | -8.32 | -8.24 | - | - | yes | yes | 53.4 | 36.6 |
| 4ty7 | 35 | 9 | 8.98 | 8.89 | 3.60 | 3.55 | -8.95 | -8.92 | - | - | no | no | 41.3 | 30.1 |
| 4u4s | 14 | 2 | 13.15 | 13.57 | 1.37 | 1.40 | -6.07 | -6.09 | - | - | yes | yes | 6.0 | 3.3 |
| 4wiv | 27 | 5 | 5.00 | 5.01 | 1.73 | 0.60 | -8.35 | -8.31 | - | - | no | no | 12.4 | 9.0 |
| 5a7b | 31 | 10 | 4.03 | 3.72 | 1.75 | 1.94 | -6.37 | -6.31 | - | - | yes | yes | 21.8 | 15.6 |
| 5aba | 23 | 5 | 9.24 | 9.13 | 1.94 | 1.90 | -6.13 | -6.08 | - | - | yes | yes | 8.4 | 6.2 |
| 5c28 | 12 | 3 | 4.72 | 4.75 | 0.82 | 1.36 | -6.39 | -6.49 | - | - | no | no | 8.6 | 5.4 |
| 5c2h | 31 | 11 | 8.93 | 7.53 | 1.88 | 1.69 | -9.54 | -9.88 | - | - | yes | yes | 22.0 | 25.8 |
| 5dwr | 32 | 6 | 1.53 | 5.55 | 1.53 | 1.51 | -10.86 | -10.96 | - | - | no | no | 20.9 | 15.5 |
| 5tmn | 32 | 17 | 10.99 | 8.51 | 2.30 | 4.49 | -8.04 | -7.67 | - | - | yes | no | 44.7 | 29.5 |

## PoseBusters failures (top pose)

- `1bzc` vina: minimum_distance_to_protein
- `1bzc` rustdock: minimum_distance_to_protein
- `1c5z` vina: minimum_distance_to_protein
- `1k1i` vina: minimum_distance_to_protein
- `1k1i` rustdock: minimum_distance_to_protein
- `1mq6` vina: minimum_distance_to_protein,volume_overlap_with_protein
- `1o0h` vina: minimum_distance_to_protein
- `1o0h` rustdock: minimum_distance_to_protein,volume_overlap_with_protein
- `1owh` vina: minimum_distance_to_protein
- `1owh` rustdock: minimum_distance_to_protein
- `1ps3` rustdock: minimum_distance_to_protein
- `1qkt` vina: minimum_distance_to_protein
- `1qkt` rustdock: minimum_distance_to_protein
- `1s38` vina: minimum_distance_to_protein,volume_overlap_with_protein
- `1s38` rustdock: minimum_distance_to_protein,volume_overlap_with_protein
- `1u1b` rustdock: minimum_distance_to_protein
- `1uto` vina: minimum_distance_to_protein
- `1uto` rustdock: minimum_distance_to_protein
- `1vso` rustdock: minimum_distance_to_protein
- `1y6r` vina: minimum_distance_to_protein
- `1y6r` rustdock: minimum_distance_to_protein
- `1ydt` vina: double_bond_flatness
- `1ydt` rustdock: double_bond_flatness
- `2cet` vina: minimum_distance_to_protein
- `2cet` rustdock: internal_energy
- `2fvd` vina: minimum_distance_to_protein
- `2fvd` rustdock: minimum_distance_to_protein
- `2fxs` vina: non-aromatic_ring_non-flatness
- `2fxs` rustdock: non-aromatic_ring_non-flatness
- `2iwx` vina: double_bond_flatness
- `2p15` vina: double_bond_flatness,minimum_distance_to_protein
- `2p15` rustdock: double_bond_flatness,minimum_distance_to_protein
- `2pog` vina: minimum_distance_to_protein
- `2pog` rustdock: minimum_distance_to_protein
- `2qbq` vina: minimum_distance_to_protein
- `2qbq` rustdock: minimum_distance_to_protein
- `2qbr` vina: minimum_distance_to_protein
- `2qbr` rustdock: minimum_distance_to_protein
- `2v7a` vina: minimum_distance_to_protein
- `2v7a` rustdock: minimum_distance_to_protein
- `2vw5` vina: double_bond_flatness,minimum_distance_to_protein
- `2vw5` rustdock: double_bond_flatness,minimum_distance_to_protein
- `2wca` vina: internal_steric_clash
- `2wca` rustdock: internal_steric_clash
- `2wer` vina: double_bond_flatness,minimum_distance_to_protein
- `2wer` rustdock: double_bond_flatness
- `2wtv` vina: minimum_distance_to_protein
- `2wtv` rustdock: minimum_distance_to_protein
- `2x00` vina: double_bond_flatness
- `2x00` rustdock: double_bond_flatness
- `2yfe` rustdock: double_bond_flatness
- `2yge` vina: double_bond_flatness,minimum_distance_to_protein
- `2yge` rustdock: double_bond_flatness,minimum_distance_to_protein
- `3arp` vina: minimum_distance_to_protein,volume_overlap_with_protein
- `3arp` rustdock: minimum_distance_to_protein,volume_overlap_with_protein
- `3arq` vina: minimum_distance_to_protein
- `3arq` rustdock: minimum_distance_to_protein
- `3b1m` vina: minimum_distance_to_protein
- `3b1m` rustdock: minimum_distance_to_protein
- `3bgz` vina: minimum_distance_to_protein
- `3bgz` rustdock: minimum_distance_to_protein
- `3coy` vina: minimum_distance_to_protein
- `3coz` vina: minimum_distance_to_protein
- `3coz` rustdock: minimum_distance_to_protein
- `3d6q` vina: minimum_distance_to_protein
- `3d6q` rustdock: minimum_distance_to_protein
- `3ebp` vina: minimum_distance_to_protein
- `3ebp` rustdock: minimum_distance_to_protein
- `3fur` vina: minimum_distance_to_protein
- `3fv1` vina: minimum_distance_to_protein
- `3fv1` rustdock: minimum_distance_to_protein
- `3fv2` vina: minimum_distance_to_protein
- `3fv2` rustdock: minimum_distance_to_protein
- `3g31` vina: minimum_distance_to_protein
- `3g31` rustdock: minimum_distance_to_protein,volume_overlap_with_protein
- `3gc5` vina: minimum_distance_to_protein,volume_overlap_with_protein
- `3gc5` rustdock: minimum_distance_to_protein,volume_overlap_with_protein
- `3gnw` vina: minimum_distance_to_protein
- `3gnw` rustdock: minimum_distance_to_protein
- `3myg` vina: minimum_distance_to_protein
- `3myg` rustdock: minimum_distance_to_protein
- `3n76` vina: double_bond_flatness
- `3n76` rustdock: double_bond_flatness
- `3oe4` vina: double_bond_flatness
- `3oe4` rustdock: double_bond_flatness
- `3oe5` rustdock: double_bond_flatness
- `3ozs` rustdock: double_bond_flatness
- `3ozt` vina: minimum_distance_to_protein
- `3p5o` vina: minimum_distance_to_protein
- `3p5o` rustdock: minimum_distance_to_protein
- `3qqs` rustdock: minimum_distance_to_protein
- `3rr4` vina: minimum_distance_to_protein,volume_overlap_with_protein
- `3rr4` rustdock: minimum_distance_to_protein,volume_overlap_with_protein
- `3ueu` vina: minimum_distance_to_protein
- `3ueu` rustdock: minimum_distance_to_protein
- `3uev` vina: minimum_distance_to_protein,volume_overlap_with_protein
- `3uev` rustdock: minimum_distance_to_protein,volume_overlap_with_protein
- `3uex` rustdock: internal_energy
- `3ui7` vina: minimum_distance_to_protein
- `3ui7` rustdock: minimum_distance_to_protein
- `3uuo` vina: minimum_distance_to_protein
- `3uuo` rustdock: minimum_distance_to_protein
- `3zt2` vina: double_bond_flatness,minimum_distance_to_protein
- `3zt2` rustdock: double_bond_flatness,minimum_distance_to_protein
- `4agp` vina: minimum_distance_to_protein
- `4agp` rustdock: minimum_distance_to_protein
- `4cig` vina: minimum_distance_to_protein
- `4cig` rustdock: minimum_distance_to_protein
- `4cra` vina: minimum_distance_to_protein
- `4cra` rustdock: minimum_distance_to_protein
- `4ddh` vina: minimum_distance_to_protein,volume_overlap_with_protein
- `4ddh` rustdock: minimum_distance_to_protein
- `4de2` rustdock: minimum_distance_to_protein
- `4dld` rustdock: minimum_distance_to_protein
- `4dli` vina: minimum_distance_to_protein
- `4e5w` vina: minimum_distance_to_protein
- `4e5w` rustdock: minimum_distance_to_protein
- `4ea2` vina: double_bond_flatness
- `4ea2` rustdock: double_bond_flatness
- `4f9w` vina: minimum_distance_to_protein
- `4f9w` rustdock: minimum_distance_to_protein
- `4gid` vina: minimum_distance_to_protein
- `4gkm` vina: minimum_distance_to_protein
- `4gkm` rustdock: minimum_distance_to_protein
- `4hge` vina: minimum_distance_to_protein
- `4j3l` vina: minimum_distance_to_protein
- `4jsz` vina: minimum_distance_to_protein
- `4jsz` rustdock: minimum_distance_to_protein
- `4k18` vina: minimum_distance_to_protein
- `4k18` rustdock: minimum_distance_to_protein
- `4m0y` vina: minimum_distance_to_protein
- `4m0y` rustdock: minimum_distance_to_protein
- `4ogj` vina: minimum_distance_to_protein,volume_overlap_with_protein
- `4ogj` rustdock: minimum_distance_to_protein,volume_overlap_with_protein
- `4ty7` vina: minimum_distance_to_protein
- `4ty7` rustdock: minimum_distance_to_protein
- `4wiv` vina: minimum_distance_to_protein
- `4wiv` rustdock: minimum_distance_to_protein
- `5c28` vina: minimum_distance_to_protein
- `5c28` rustdock: minimum_distance_to_protein,volume_overlap_with_protein
- `5dwr` vina: minimum_distance_to_protein
- `5dwr` rustdock: minimum_distance_to_protein
- `5tmn` rustdock: minimum_distance_to_protein
