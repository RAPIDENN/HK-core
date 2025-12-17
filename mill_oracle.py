#!/usr/bin/env python3
import argparse
import json
import subprocess
import time
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser("mill_oracle")
    p.add_argument("--url", default="http://127.0.0.1:8080/mill/refine", help="/mill/refine URL")
    p.add_argument("--token", default="devtoken", help="Bearer token")
    p.add_argument("--ls", default="8,16,32,64", help="Lattice sizes, e.g. 8,16,32,64")
    p.add_argument("--beta", type=float, default=2.0, help="Gauge coupling beta")
    p.add_argument("--seed", type=int, default=123, help="Random seed")
    p.add_argument("--n-thermal-sweeps", type=int, default=200, help="Thermalization sweeps")
    p.add_argument("--n-sweeps", type=int, default=1000, help="Initial sweeps")
    p.add_argument("--measure-every", type=int, default=10, help="Measurement period")
    p.add_argument("--step-size", type=float, default=0.3, help="Metropolis step size")
    p.add_argument("--max-rounds", type=int, default=3, help="Max retry rounds")
    p.add_argument("--sleep", type=float, default=0.5, help="Sleep seconds between rounds")
    p.add_argument("--output", required=True, help="Output JSON file path")
    p.add_argument("--quiet", action="store_true", help="Minimal stdout (progress only)")
    return p.parse_args()


def call_refine(url, token, payload):
    cmd = [
        "curl",
        "-sS",
        "-X",
        "POST",
        url,
        "-H",
        f"Authorization: Bearer {token}",
        "-H",
        "Content-Type: application/json",
        "-d",
        json.dumps(payload),
    ]
    out = subprocess.check_output(cmd)
    return json.loads(out)


def main():
    args = parse_args()
    ls = [int(x) for x in args.ls.split(",") if x.strip()]

    history = []
    n_sweeps = args.n_sweeps
    stop_reason = "max_rounds"

    for r in range(1, args.max_rounds + 1):
        payload = {
            "ls": ls,
            "beta": args.beta,
            "n_thermal_sweeps": args.n_thermal_sweeps,
            "n_sweeps": n_sweeps,
            "measure_every": args.measure_every,
            "step_size": args.step_size,
            "seed": args.seed,
        }

        res = call_refine(args.url, args.token, payload)
        fv = res["result"]["final_verdict"]
        status = fv["status"]

        history.append({
            "round": r,
            "n_sweeps": n_sweeps,
            "final_verdict": fv,
        })

        print(f"[round {r}] n_sweeps={n_sweeps} status={status}")

        if status != "inconclusive":
            stop_reason = "decision"
            break

        if r < args.max_rounds:
            n_sweeps *= 2
            time.sleep(args.sleep)

    summary = {
        "config": {
            "url": args.url,
            "ls": ls,
            "beta": args.beta,
            "seed": args.seed,
            "n_thermal_sweeps": args.n_thermal_sweeps,
            "measure_every": args.measure_every,
            "step_size": args.step_size,
            "initial_n_sweeps": args.n_sweeps,
            "max_rounds": args.max_rounds,
        },
        "rounds": history,
        "stopped_because": stop_reason,
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(summary, indent=2))

    if not args.quiet:
        print("\nSUMMARY written to", out_path)


if __name__ == "__main__":
    main()
