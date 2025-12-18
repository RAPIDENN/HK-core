# Architecture (su2-3d)

This branch keeps the same HTTP API and JSON schema as `master` and swaps the 3D kernel to **SU(2) in 3D** while keeping the analysis/oracle pipeline intact.

## Data Flow

1. `POST /mill/refine` (API)
2. Per-`L` run loop:
   - initialize field (IC)
   - thermal sweeps
   - measurement sweeps
3. Observables + analysis:
   - plaquette stats
   - invariance scaling + reflection positivity telemetry
   - effective mass estimator + plateau detection
   - final verdict synthesis (unchanged)

## Opt-in runtime controls (env vars)

These do not change the request schema and are meant for controlled experiments.

- **Initial conditions**: `MILL_IC`
  - `cold` (default): identity links
  - `hot`: Haar-random links
  - `smooth`: hot start + pre-smoothing (APE) before thermalization

- **Update dynamics**: `MILL_UPDATE`
  - unset (default): Metropolis random-link sweep
  - `hb` / `heatbath`: SU(2) heatbath sweep
  - `hb_or` / `heatbath_overrelax`: heatbath sweep + overrelaxation sweep

## Tracked results

Canonical outputs for this branch live under `results/`.
- `results/README.md`

