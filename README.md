# HK-Core: Autonomous Lattice Field Theory Engine

Autonomous lattice gauge simulation engine in Rust. Validated baseline: 2D U(1). Current experimental kernel: 4D SU(3) pure gauge with Cabibbo-Marinari SU(2)-subgroup heatbath updates.

> **Status:** 2D U(1) baseline validated at `L=64`; 4D SU(3) is experimental and has no final mass-gap claim.
> **License:** GNU AGPLv3

## Overview
HK-Core is a high-performance, autonomous simulation engine written in **Rust**. It keeps the validated **2D U(1) Lattice Gauge Theory** kernel as a regression baseline and now includes an experimental **4D SU(3) pure gauge** kernel for quenched Yang-Mills studies.

Unlike traditional academic codes that output raw data for manual analysis, HK-Core acts as an **autonomous scientific instrument**. It incorporates an internal decision engine that performs statistical blocking, self-calibrates via **Jackknife resampling**, and issues binary threshold verdicts without human intervention.

## Key Capabilities
1. **Autonomous Decision Making:** Implements an **IR Dominance** strategy (`verdict_mode="ir_lmax"`), prioritizing the infrared signal at the largest requested lattice size.
2. **SU(3) Pure Gauge Kernel:** Runs 4D SU(3) links without dynamical fermions. The normal SU(3) update path is `MILL_UPDATE=hb_or`, a Cabibbo-Marinari heatbath sweep followed by overrelaxation; `MILL_UPDATE=hb` is available for pure heatbath, and `metropolis`/`metro_or` are explicit comparison modes.
3. **Rigorous Statistics:** Native implementation of **APE Smearing**, effective-mass plateau extraction, and **Jackknife Blocking** to isolate signals from Monte Carlo noise.
4. **Physical Integrity:** Runtime diagnostics for **Reflection Positivity** and **Translation Invariance** scaling.
5. **Reproducibility:** Deterministic baseline outcomes verified across independent seeds (`777`, `111`, `1234`), plus tracked small regression runs under `results/`.

## Validation
For the validated configuration (`beta=2.0`, `step_size=0.06`, `n_thermal_sweeps=2000`, `n_sweeps=60000`, `measure_every=10`, `ls=[8,16,32,64]`, `plateau_mode="stat"`, `verdict_mode="ir_lmax"`), seeds `777`, `111`, `1234` reproduce:
- `final_verdict.ir_lmax.per_m0["m0=0.1"] == "compatible"`

Interpretation: `"compatible"` for `m0=0.1` means `mean - k_sigma*std > 0.1` with `k_sigma=2.0`, `plateau_width>=6`, and (in `plateau_mode="stat"`) a `chi2_ok=true` safeguard.

See **Releases** for canonical JSON datasets + checksums (outputs are not committed under `out/`).

## SU(3) status
The SU(3) path is intentionally conservative:
- The theory is pure gauge: there are no quarks or fermion determinants in the kernel.
- The default update mode is Cabibbo-Marinari heatbath plus overrelaxation (`MILL_UPDATE=hb_or`). Plain Metropolis must be requested explicitly and should be treated as a diagnostic/comparison mode, not the default production path.
- Plateau extractors reject negative effective-mass windows. A negative window is treated as no physical mass plateau, not as a small-gap signal.
- Recent beta=5.5 heatbath sweeps produced promising short-window mass estimates, but longer and larger-volume checks did not sustain a robust plateau. Do not cite those runs as a closed mass-gap result.

## Tracked results
Small canonical outputs are also tracked under `results/` for quick inspection:
- 2D baseline (seeds 777/111/1234): `results/2d_u1/beta2_seed777_s60000_stat_ir_lmax.json`
- 3D references (U(1) dim-check, SU(2) IC/IR instrumentation): `results/README.md`
- 4D SU(3) regression references: `results/README.md`

## Architecture
See `ARCHITECTURE.md` for the system map (pipelines, contracts, branch strategy).

## Run locally
```bash
export AUTH_TOKEN=devtoken
export RUST_LOG=info
export MILL_UPDATE=hb_or
cargo run --release
```

`MILL_UPDATE` may be set to `hb`, `hb_or`, `metropolis`, `metro_or`, or `or`. If unset, SU(3) defaults to `hb_or`.

## SU(3) 4D smoke
```bash
curl -sS -X POST "http://127.0.0.1:8080/mill/refine" \
  -H "Authorization: Bearer devtoken" \
  -H "Content-Type: application/json" \
  -d '{
    "ls": [8,12],
    "beta": 5.5,
    "n_thermal_sweeps": 120,
    "n_sweeps": 240,
    "measure_every": 10,
    "step_size": 0.2,
    "seed": 944,
    "verdict_mode": "ir_lmax",
    "plateau_mode": "stat"
  }' \
| python3 -c 'import json,sys; d=json.load(sys.stdin)["result"]["final_verdict"]; ir=d.get("ir_lmax",{}); print(d["status"], ir.get("channel"), ir.get("plateau_width"), ir.get("m_eff_mean"), ir.get("m_eff_std"))'
```

## IR-Lmax check (m0=0.1)
```bash
curl -sS -X POST "http://127.0.0.1:8080/mill/refine" \
  -H "Authorization: Bearer devtoken" \
  -H "Content-Type: application/json" \
  -d '{
    "ls": [8,16,32,64],
    "beta": 2.0,
    "n_thermal_sweeps": 2000,
    "n_sweeps": 60000,
    "measure_every": 10,
    "step_size": 0.06,
    "seed": 777,
    "verdict_mode": "ir_lmax",
    "plateau_mode": "stat"
  }' \
| python3 -c 'import json,sys; d=json.load(sys.stdin); v=d["result"]["final_verdict"]; ir=v.get("ir_lmax",{}); per=ir.get("per_m0",{}); print("IR_LMAX m0=0.1:", per.get("m0=0.1")); print("L:", ir.get("l"), "channel:", ir.get("channel"), "width:", ir.get("plateau_width"), "mean:", ir.get("m_eff_mean"), "std:", ir.get("m_eff_std"), "chi2_ok:", ir.get("chi2_ok"))'
```

## Repo hygiene
- Build artifacts live in `target/` and are ignored via `.gitignore`.
- Experiment outputs live in `out/` and are ignored; publish “canonical” outputs as GitHub release assets instead.

## License
Licensed under the GNU AGPLv3. See `LICENSE`.
