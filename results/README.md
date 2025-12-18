# Results (tracked)

This branch (`su2-3d`) keeps a small set of **tracked, canonical JSON outputs** under `results/` for auditability.

Notes:
- These files are produced by local runs of `POST /mill/refine` and are meant to be small and stable.
- Large, exploratory batches should stay in `out/` (ignored by git).

## 3D SU(2) — A/B dynamics check (Metropolis vs HB+OR)

Payload (shared across runs):
- `results/3d_su2/ab_payload_beta1.6_ls16-32_th2000_sw8000_me10_stat_ir_lmax.json`

Outputs:
- Baseline (Metropolis, `MILL_UPDATE` unset): `results/3d_su2/ab_seed777_metropolis_beta1.6_ls16-32_th2000_sw8000_me10_stat_ir_lmax.json`
- HB+OR (Heatbath+Overrelax, `MILL_UPDATE=hb_or`): `results/3d_su2/ab_seed777_hb_or_beta1.6_ls16-32_th2000_sw8000_me10_stat_ir_lmax.json`
- HB+OR seeds (stability check):
  - `results/3d_su2/ab_seed111_hb_or_beta1.6_ls16-32_th2000_sw8000_me10_stat_ir_lmax.json`
  - `results/3d_su2/ab_seed1234_hb_or_beta1.6_ls16-32_th2000_sw8000_me10_stat_ir_lmax.json`

