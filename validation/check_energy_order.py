#!/usr/bin/env python3
"""Compare energy error and step-size scaling on identical initial states.

Runs all three implementations for the same physical duration. OpenMP uses
one thread as a serial control and --threads for the parallel comparison.
Raw jobs are temporary; CSV summaries, traces, figures and a short report
are saved together. Requires NumPy, Matplotlib and the C++ build dependencies.
"""

import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
VARIANTS = ("MDNormal", "MDMultipleTimeSteps", "MDMultipleTimeStepsOpenMP")
LABELS = dict(zip(VARIANTS, ("Uniform Verlet", "Serial MTS", "OpenMP MTS")))


def run(command, **kwargs):
    result = subprocess.run(list(map(str, command)), capture_output=True, text=True,
                            timeout=180, **kwargs)
    if result.returncode:
        raise RuntimeError(f"Command failed: {command}\n{result.stdout}\n{result.stderr}")
    return result.stdout


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def statistics(path, dt, variant, threads):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    energy = np.array([float(row["Etot"]) for row in rows])
    times = np.array([int(row["t"]) * dt for row in rows])
    if len(energy) < 5 or not np.isfinite(energy).all():
        raise ValueError(f"Invalid energy trace: {path}")
    # Equal-duration sampling: use the second half of [0, duration), omitting
    # the final endpoint for the fluctuation fit. Error extrema use all rows.
    sample = energy[:-1]
    tail = sample[len(sample)//2:]
    tail_times = times[:-1][len(sample)//2:]
    deviation = float(np.std(tail, ddof=0))
    summary = {
        "variant": variant, "threads": threads, "dt": dt,
        "duration": float(times[-1]), "rows": len(rows),
        "E0": float(energy[0]), "K0": float(rows[0]["K"]), "U0": float(rows[0]["U"]),
        "std_energy": deviation, "variance_energy": deviation**2,
        "max_abs_energy_error": float(np.max(np.abs(energy-energy[0]))),
        "final_energy_error": float(energy[-1]-energy[0]),
        "tail_linear_drift_per_time": float(np.polyfit(tail_times, tail, 1)[0]),
    }
    traces = [{"variant": variant, "threads": threads, "dt": dt,
               "time": float(t), "energy_error": float(e-energy[0])}
              for t, e in zip(times, energy)]
    return summary, traces


def fit(rows, variant, threads):
    selected = sorted((r for r in rows if r["variant"] == variant and r["threads"] == threads),
                      key=lambda r: r["dt"])
    steps = np.array([r["dt"] for r in selected])
    deviations = np.array([r["std_energy"] for r in selected])
    if np.any(deviations <= 0):
        raise ValueError("Need nonzero energy fluctuations to fit an order")
    slope, intercept = np.polyfit(np.log(steps), np.log(deviations), 1)
    observed = np.log(deviations)
    residual = observed - (slope*np.log(steps)+intercept)
    return {"variant": variant, "threads": threads, "std_slope": float(slope),
            "variance_slope": float(2*slope),
            "r_squared_log_std": float(1-np.sum(residual**2)/np.sum((observed-observed.mean())**2))}


def label(variant, threads):
    suffix = f" ({threads} thread{'s' if threads != 1 else ''})"
    return LABELS[variant] + (suffix if variant.endswith("OpenMP") else "")


def figures(rows, traces, fits, output, ratio, parallel_threads):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5))
    for result in fits:
        variant, threads = result["variant"], result["threads"]
        # The one-thread control remains in the tables/CSVs. Plot the three
        # implementations without a fourth curve hiding the serial MTS line.
        if variant.endswith("OpenMP") and threads != parallel_threads:
            continue
        selected = sorted((r for r in rows if r["variant"] == variant and r["threads"] == threads),
                          key=lambda r: r["dt"])
        ax.loglog([r["dt"] for r in selected], [r["std_energy"] for r in selected], "o-",
                  label=f"{label(variant, threads)}: slope {result['std_slope']:.2f}")
    serial = sorted((r for r in rows if r["variant"] == VARIANTS[1]), key=lambda r:r["dt"])
    ax.loglog([r["dt"] for r in serial],
              [serial[0]["std_energy"]*(r["dt"]/serial[0]["dt"])**2 for r in serial],
              "k--", alpha=0.5, label="Quadratic reference (slope 2)")
    ax.set_xlabel(f"Macro step h (small step h/{ratio})")
    ax.set_ylabel("Energy standard deviation, second half of run")
    ax.set_title("Same initial state and physical duration")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(output/"energy_order.png", dpi=170)
    plt.close(fig)

    steps = sorted({r["dt"] for r in rows})
    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    for ax, dt in zip(axes, (steps[0], steps[-1])):
        for variant in VARIANTS:
            threads = parallel_threads if variant.endswith("OpenMP") else 1
            selected = [r for r in traces if r["variant"] == variant and r["threads"] == threads and r["dt"] == dt]
            ax.plot([r["time"] for r in selected], [r["energy_error"] for r in selected],
                    label=label(variant, threads), linewidth=0.8)
        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_title(f"Macro step h = {dt:g}")
        ax.set_ylabel("E(t) - E(0)")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    axes[-1].set_xlabel("Physical time")
    fig.tight_layout()
    fig.savefig(output/"energy_error.png", dpi=170)
    plt.close(fig)


def report(output, metadata, rows, fits):
    algorithm_link = os.path.relpath(ROOT/"docs/ALGORITHM.md", output.resolve())
    lines = ["# Energy error and integration order", "",
             "All implementations start from the same particle CSV and use the same",
             "configuration, macro step and physical duration. The build is `releaseO3`.",
             "Uniform Verlet updates every particle with `h/Ratio`; MTS updates the",
             "heavy group with `h` and holds the other group fixed in each stage.", "",
             f"Duration: {metadata['duration']:g}. Step sizes: " +
             ", ".join(f"`{h:g}`" for h in metadata["step_sizes"]) + ".", "",
             "## Energy fluctuations", "",
             "A slope near 2 means the standard deviation decreases approximately",
             "quadratically with the step. Its variance then decreases with the fourth",
             "power. The fit uses the second half of the run, excluding the final endpoint.", "",
             "| Implementation | Standard-deviation slope | Variance slope |",
             "| --- | ---: | ---: |"]
    for row in fits:
        lines.append(f"| {label(row['variant'], row['threads'])} | {row['std_slope']:.3f} | {row['variance_slope']:.3f} |")
    lines += ["", "![Energy fluctuations versus step](energy_order.png)", "",
              "## Error relative to the initial energy", "",
              "![Physical energy error](energy_error.png)", "",
              "| Macro step | Implementation | Maximum abs(E - E0) | Final E - E0 |",
              "| --- | --- | ---: | ---: |"]
    for row in rows:
        lines.append(f"| {row['dt']:g} | {label(row['variant'], row['threads'])} | {row['max_abs_energy_error']:.6g} | {row['final_energy_error']:+.6g} |")
    lines += ["", "Both methods can be second order with different error amplitudes.",
              "The uniform baseline resolves heavy motion with a smaller step too.",
              "MTS also introduces error by freezing the inactive group during a stage.",
              "The relative error amplitudes depend on the state, mass ratio and step;",
              "there is no universal MTS/uniform error factor.", "",
              "One-thread OpenMP is the control for the serial MTS update. Multiple",
              "threads can produce different long trajectories through force-summation",
              "roundoff. Finite-run maxima are not a systematic parallel accuracy ratio.", "",
              "These fits and trajectories do not prove a long-time energy bound.",
              "The force is discontinuous at the cutoff, which limits smooth-force",
              f"convergence arguments. See [the algorithm guide]({algorithm_link}).", "",
              "## Reproduce", "", "From the repository root:", "", "```bash",
              "python3 validation/check_energy_order.py", "```", "",
              "Use `--eigen-include /path/to/eigen` for a custom installation. Options",
              "`--dt`, `--duration`, `--threads` and `--output-dir` select the experiment.",
              "Raw jobs are temporary. [results.json](results.json) records the input",
              "hash, configuration and fits; [energy_statistics.csv](energy_statistics.csv)",
              "and [energy_traces.csv](energy_traces.csv) contain the numerical results.", ""]
    (output/"README.md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT/"validation/results/energy_order")
    parser.add_argument("--eigen-include", default=os.environ.get("EIGEN_INC", "/usr/include/eigen3"))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--dt", nargs="+", type=float, default=[0.04, 0.05, 0.08, 0.1])
    parser.add_argument("--duration", type=float, default=50)
    args = parser.parse_args()
    if args.threads <= 0: parser.error("Thread count must be positive")
    steps = sorted(set(args.dt))
    if len(steps) < 2 or any(not np.isfinite(h) or h <= 0 for h in steps):
        parser.error("Use at least two distinct positive finite step sizes")
    if not np.isfinite(args.duration) or args.duration <= 0:
        parser.error("Duration must be positive and finite")
    if any(round(args.duration/h) < 4 or not np.isclose(round(args.duration/h)*h, args.duration) for h in steps):
        parser.error("Duration must be a whole number of at least four steps for every dt")
    if not (Path(args.eigen_include)/"Eigen/Core").is_file(): parser.error("Eigen/Core not found")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location("generate_ic", ROOT/VARIANTS[0]/"jobs/sim1/generateIC.py")
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)
    cfg = generator.load_config(ROOT/VARIANTS[0]/"jobs/sim1/config.txt")
    ratio = max(1, int(np.sqrt(max(cfg["mass1"], cfg["mass2"])/min(cfg["mass1"], cfg["mass2"]))))
    initial = ROOT/VARIANTS[0]/"jobs/sim1/initial_conf.csv"
    cases = [(VARIANTS[0], 1), (VARIANTS[1], 1), (VARIANTS[2], 1)]
    if args.threads != 1: cases.append((VARIANTS[2], args.threads))
    rows, traces = [], []
    for variant in VARIANTS:
        run(["make", "-C", ROOT/variant, "-j2", "BUILD_MODE=releaseO3", "EIGEN_INC="+args.eigen_include])
    with tempfile.TemporaryDirectory(prefix="md-energy-") as temp:
        for variant, threads in cases:
            for h in steps:
                job = Path(temp)/f"{variant}_{threads}_{h}"
                job.mkdir()
                config = dict(cfg, dt=h, totsteps=round(args.duration/h))
                (job/"config.txt").write_text("".join(f"{key} {value}\n" for key, value in config.items()))
                shutil.copyfile(initial, job/"initial_conf.csv")
                run([ROOT/variant/"bin/MolecularDynamics_O3", job],
                    env=dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE"))
                row, trace = statistics(job/"energy.csv", h, variant, threads)
                rows.append(row)
                traces.extend(trace)
                print(f"{label(variant, threads)}, h={h:g}: std(E)={row['std_energy']:.9g}", flush=True)
    fits = [fit(rows, variant, threads) for variant, threads in cases]
    metadata = {"duration": args.duration, "step_sizes": steps, "tail_fraction": 0.5,
                "mass_ratio_substeps": ratio, "build_mode": "releaseO3",
                "configuration": {key:value for key,value in cfg.items() if key not in ("dt", "totsteps")},
                "initial_conditions": "MDNormal/jobs/sim1/initial_conf.csv",
                "initial_conditions_sha256": hashlib.sha256(initial.read_bytes()).hexdigest(),
                "compiler": run([*shlex.split(os.environ.get("CXX", "g++")), "--version"]).splitlines()[0]}
    (args.output_dir/"results.json").write_text(json.dumps({"metadata":metadata, "fits":fits,
                                                        "energy_statistics":rows}, indent=2, allow_nan=False)+"\n")
    write_csv(args.output_dir/"energy_statistics.csv", rows)
    write_csv(args.output_dir/"energy_traces.csv", traces)
    figures(rows, traces, fits, args.output_dir, ratio, args.threads)
    report(args.output_dir, metadata, rows, fits)
    print(json.dumps(fits, indent=2))


if __name__ == "__main__":
    main()
