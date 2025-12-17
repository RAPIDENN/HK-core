# SU(2) 3D experimental branch

This branch (`su2-3d`) is an **experimental research line**.

## What it is
- Replaces the 3D U(1) kernel with an **SU(2) 3D** kernel (plaquettes, Metropolis updates, APE smearing with explicit SU(2) projection).
- Keeps **API, JSON schema, oracle, classifiers, and knobs unchanged**.
- Intended for **instrument validation and observability-limit characterization**, not for physics claims.

## Instrumental IC (coupling test)
This branch adds an opt-in initial-condition switch via environment variable:
- `MILL_IC=cold` (default): cold/identity start
- `MILL_IC=hot`: Haar-random hot start
- `MILL_IC=smooth`: hot start + pre-smoothing (APE) before thermalization

These modes are for **instrumental coupling tests** only. Results obtained under non-default IC must be reported as such.

## Run locally
```bash
export AUTH_TOKEN=devtoken
export RUST_LOG=info
# optional:
# export MILL_IC=hot
# export MILL_IC=smooth
cargo run --release
```

## Notes
- `master` remains the canonical 2D U(1) validated release line.
- Outputs in `out/` are intentionally not committed.
