# HK-Core: Autonomous Lattice Field Theory Engine

> **Status:** VALIDATED (Autonomous mass-gap-like threshold check confirmed at `L=64` (2D))
> **License:** GNU AGPLv3

## ⚡ Overview
HK-Core is a high-performance, autonomous simulation engine written in **Rust**. It implements a **2D U(1) Lattice Gauge Theory** kernel as a proof-of-concept to demonstrate an industrial-grade architecture for automated scientific discovery.

Unlike traditional academic codes that output raw data for manual analysis, HK-Core acts as an **autonomous scientific instrument**. It incorporates an internal decision engine that performs statistical blocking, self-calibrates via **Jackknife resampling**, and issues binary threshold verdicts without human intervention.

## 🚀 Key Capabilities
1. **Autonomous Decision Making:** Implements an **IR Dominance** strategy (`verdict_mode="ir_lmax"`), prioritizing the infrared signal at `L=64` (2D lattice) to filter out finite-volume noise.
2. **Rigorous Statistics:** Native implementation of **APE Smearing** and **Jackknife Blocking** to isolate signals from Monte Carlo noise.
3. **Physical Integrity:** Runtime diagnostics for **Reflection Positivity** and **Translation Invariance** scaling.
4. **Reproducibility:** Deterministic outcomes verified across independent seeds (`777`, `111`, `1234`).

## 📊 Validation
For the validated configuration (`beta=2.0`, `step_size=0.06`, `n_thermal_sweeps=2000`, `n_sweeps=60000`, `measure_every=10`, `ls=[8,16,32,64]`, `plateau_mode="stat"`, `verdict_mode="ir_lmax"`), seeds `777`, `111`, `1234` reproduce:
- `final_verdict.ir_lmax.per_m0["m0=0.1"] == "compatible"`

Interpretation: `"compatible"` for `m0=0.1` means `mean - k_sigma*std > 0.1` with `k_sigma=2.0`, `plateau_width>=6`, and (in `plateau_mode="stat"`) a `chi2_ok=true` safeguard.

See **Releases** for canonical JSON datasets + checksums (outputs are not committed under `out/`).

## 🛠️ Run locally
```bash
export AUTH_TOKEN=devtoken
export RUST_LOG=info
cargo run --release
```

## 🧪 IR-Lmax check (m0=0.1)
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
| python3 - <<'PY'
import sys, json
d=json.load(sys.stdin)
v=d["result"]["final_verdict"]
ir=v.get("ir_lmax",{})
per=ir.get("per_m0",{})
print("IR_LMAX m0=0.1:", per.get("m0=0.1"))
print("L:", ir.get("l"), "channel:", ir.get("channel"), "width:", ir.get("plateau_width"),
      "mean:", ir.get("m_eff_mean"), "std:", ir.get("m_eff_std"), "chi2_ok:", ir.get("chi2_ok"))
PY
```

## Repo hygiene
- Build artifacts live in `target/` and are ignored via `.gitignore`.
- Experiment outputs live in `out/` and are ignored; publish “canonical” outputs as GitHub release assets instead.

## License
Licensed under the GNU AGPLv3. See `LICENSE`.
