# Results (tracked)

This folder contains a small set of **tracked** JSON outputs used as a reproducible reference.

Notes:
- Large/iterative outputs still belong in `out/` (gitignored).
- These files are intentionally small and few.

## 2D U(1) baseline (master)
- `results/2d_u1/beta2_seed777_s60000_stat_ir_lmax.json`
- `results/2d_u1/beta2_seed111_s60000_stat_ir_lmax.json`
- `results/2d_u1/beta2_seed1234_s60000_stat_ir_lmax.json`

## 3D U(1) (transition-legacy-3d)
- `results/3d_u1/demo_3d_u1_run_L16.json` (dim check via `/mill/run`)

## 3D SU(2) (su2-3d)
These are **instrumental** reference runs (not physics claims), included to document the current max IR reached and IC coupling tests:
- `results/3d_su2/ic_test_A_hot_seed777_ls16-32_b1p6_s0p06_th2000_sw8000_me10_stat_ir_lmax.json` (best observed `ir_lmax.plateau_width` so far)
- `results/3d_su2/ord021d_su2_seed777_ls16-32_b1p6_s0p06_th2000_sw8000_me10_stat_ir_lmax.json` (baseline thermal IC)

## 4D SU(3) pure gauge (master)
These are small regression references for the current SU(3) kernel. They are **not** final mass-gap claims.
- `results/4d_su3/beta5p5_seed944_heatbath_positive_plateau_regression.json` verifies that negative/invalid effective-mass windows do not count as physical plateaus after the positive-mass plateau filter. The run is `MILL_UPDATE=heatbath`, `ls=[8,12]`, `beta=5.5`, `thermal=120`, `sweeps=240`, `seed=944`, and remains `final_status=inconclusive`.

## Index
- `results/RESULTS_INDEX.json` contains a machine-readable summary of the tracked runs.
