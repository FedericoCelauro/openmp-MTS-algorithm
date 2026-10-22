#!/usr/bin/env python3
"""Compare real integration kernels and complete applications on identical jobs.

Builds all variants with releaseO3. Each case has one generated initial state
shared by every implementation and thread count. Runs are sequential, with a
warmup, randomized repetition order, and medians/ranges. Jobs are temporary.
Speedups compare equal macro steps and physical duration, not equal error.
Optional NumPy/Matplotlib dependencies are in requirements-analysis.txt.
"""

import argparse
import contextlib
import csv
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import platform
import random
import shlex
import shutil
import statistics
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
VARIANTS = ("MDNormal", "MDMultipleTimeSteps", "MDMultipleTimeStepsOpenMP")
CASES = ((4, 0.1, 200), (8, 0.01, 30), (8, 0.1, 30), (8, 0.5, 30))


def run(command, **kwargs):
    result = subprocess.run(list(map(str, command)), capture_output=True, text=True,
                            timeout=180, **kwargs)
    if result.returncode:
        raise RuntimeError(f"Command failed: {command}\n{result.stdout}\n{result.stderr}")
    return result


def load_generator():
    spec = importlib.util.spec_from_file_location("generate_ic", ROOT/VARIANTS[0]/"jobs/sim1/generateIC.py")
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)
    return generator


def write_config(job, cfg):
    (job/"config.txt").write_text("\n".join(f"{key} {value}" for key, value in cfg.items())+"\n")


def pair_work(count, fraction, ratio=10):
    """Count pair visits per macro step, including rejected cutoff pairs.

    Applicable to the even-ratio scheme used in this benchmark. This is an
    operation-count prediction, not a model for atomic/thread or memory costs.
    """
    light = int(count*fraction)
    heavy = count-light
    choose = lambda n: n*(n-1)//2
    uniform = ratio*choose(count)
    mts = ratio*(choose(light)+light*heavy) + choose(heavy)+light*heavy
    return {"light_count": light, "heavy_count": heavy,
            "uniform_pair_visits": uniform, "mts_pair_visits": mts,
            "pair_visit_reduction": uniform/mts}


def build(eigen, temporary):
    binaries = {}
    for variant in VARIANTS:
        project = ROOT/variant
        run(["make", "-C", project, "-j2", "BUILD_MODE=releaseO3", "EIGEN_INC="+eigen])
        flags = ["-std=c++17", "-O3", "-I"+str(project/"src/include"), "-I"+eigen]
        if variant != "MDNormal": flags.append("-DMTS")
        if variant.endswith("OpenMP"): flags.append("-fopenmp")
        objects = sorted(p for p in (project/"build/releaseO3").rglob("*.o") if p.name != "MolecularDynamics.o")
        binary = temporary/(variant+"_kernel")
        run([*shlex.split(os.environ.get("CXX", "g++")), *flags, ROOT/"validation/support/benchmark_integrators.cpp",
             *objects, "-o", binary, "-lstdc++fs"])
        binaries[variant] = binary
    return binaries


def environment(threads):
    return dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE",
                OMP_PROC_BIND="close", OMP_PLACES="cores", OMP_WAIT_POLICY="ACTIVE")


def measure(mode, variant, threads, job, binaries):
    env = environment(threads)
    if mode == "kernel":
        data = json.loads(run([binaries[variant], job], env=env).stdout)
        if data["threads"] != threads:
            raise RuntimeError(f"Requested {threads} threads, got {data['threads']}")
        return data["seconds"]
    start = time.perf_counter()
    run([ROOT/variant/"bin/MolecularDynamics_O3", job], env=env)
    elapsed = time.perf_counter()-start
    # Excluded from elapsed: validate successful, finite physical output.
    with (job/"energy.csv").open() as stream:
        energy = list(csv.DictReader(stream))
    if not all(math.isfinite(float(value)) for row in energy for value in row.values()):
        raise RuntimeError("Non-finite application energy")
    return elapsed


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_speed(rows, output, configurations):
    import numpy as np
    import matplotlib.pyplot as plt
    cases = list(dict.fromkeys(row["case"] for row in rows))
    labels = [f"N={row['particles']}\nlight={row['light_fraction']:.0%}"
              for case in cases for row in rows if row["case"]==case and row["mode"]=="kernel" and row["variant"]==VARIANTS[0]]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    width = 0.8/len(configurations)
    for ax, mode in zip(axes, ("kernel", "application")):
        for index, (variant, threads) in enumerate(configurations):
            selected = [next(row for row in rows if row["case"]==case and row["mode"]==mode
                             and row["variant"]==variant and row["threads"]==threads) for case in cases]
            values = [row["speedup_vs_uniform"] for row in selected]
            label = "Uniform" if variant==VARIANTS[0] else "Serial MTS" if variant==VARIANTS[1] else f"OpenMP {threads}T"
            lower, upper = [], []
            for row in selected:
                baseline = next(r["median_seconds"] for r in rows if r["case"]==row["case"]
                                and r["mode"]==mode and r["variant"]==VARIANTS[0])
                lower.append(row["speedup_vs_uniform"]-baseline/row["max_seconds"])
                upper.append(baseline/row["min_seconds"]-row["speedup_vs_uniform"])
            ax.bar(np.arange(len(cases))+(index-(len(configurations)-1)/2)*width, values, width,
                   yerr=np.array([lower, upper]), capsize=2, label=label,
                   error_kw={"linewidth": 0.8})
        ax.set_xticks(np.arange(len(cases)))
        ax.set_xticklabels(labels)
        ax.axhline(1, color="black", linewidth=0.8)
        ax.set_title("Integration only" if mode=="kernel" else "Complete application, including output")
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("Speedup relative to uniform integration (larger is faster)")
    axes[-1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output/"speed_comparison.png", dpi=170)
    plt.close(fig)


def report(output, metadata, summaries, cases):
    """Keep the readable report synchronized with the generated timing data."""
    algorithm_link = os.path.relpath(ROOT/"docs/ALGORITHM.md", output.resolve())
    energy_link = os.path.relpath(ROOT/"validation/results/energy_order/README.md", output.resolve())
    names = dict(zip(VARIANTS, ("Uniform Verlet", "Serial MTS", "OpenMP MTS")))
    representative = next(case for case in cases if case["particles"] == 512 and case["light_count"] == 51)
    selected = [row for row in summaries if row["case"] == representative["case"]]
    configurations = [(row["variant"], row["threads"]) for row in selected if row["mode"] == "application"]
    lines = ["# Speed comparison", "",
             "All implementations use identical initial conditions within each case,",
             "the same macro step and physical duration, and a `releaseO3` build.",
             f"Each configuration has one warmup and {metadata['repeats']} timed repetitions.",
             "Runs are sequential in randomized order. Reported times are medians;",
             "the CSV files retain the individual samples and timing ranges.", "",
             f"CPU: {metadata['cpu_model']}. Compiler: {metadata['compiler']}.",
             "OpenMP uses fixed team sizes, core affinity and active waiting; full",
             "settings are recorded in [results.json](results.json).", "",
             "## Representative case", "",
             "512 particles, 51 light particles, masses 1 and 100, macro step 0.01,",
             "30 macro steps (physical duration 0.3):", "",
             "| Implementation | Integration seconds | Integration speedup | Complete seconds | Complete speedup |",
             "| --- | ---: | ---: | ---: | ---: |"]
    for variant, threads in configurations:
        kernel = next(r for r in selected if r["variant"] == variant and r["threads"] == threads and r["mode"] == "kernel")
        application = next(r for r in selected if r["variant"] == variant and r["threads"] == threads and r["mode"] == "application")
        name = names[variant] + (f" ({threads} threads)" if variant.endswith("OpenMP") else "")
        lines.append(f"| {name} | {kernel['median_seconds']:.4f} | {kernel['speedup_vs_uniform']:.2f}x | {application['median_seconds']:.4f} | {application['speedup_vs_uniform']:.2f}x |")
    lines += ["", "Integration timing excludes initialization, energy logging and snapshots.",
              "Complete timing includes process startup and output with warm filesystem",
              "caches. Worker startup is excluded only from the integration timer.", "",
              "## Size and mixture dependence", "",
              "![Integration and complete-application speedups](speed_comparison.png)", "",
              "MTS recalculates heavy–heavy forces once per macro step, rather than",
              "at every small step. Its benefit grows when few particles are light.",
              "OpenMP adds synchronization and atomic-update overhead, so small cases",
              "can benefit less from parallelism. Candidate-pair counts for each case",
              "are included in the JSON metadata.", "",
              "The 64-particle case uses 200 macro steps; 512-particle cases use 30.",
              "Compare implementations within a case; absolute times across cases",
              "with different durations are not a size-scaling experiment.", "",
              "These are fixed-step speedups, not matched-accuracy gains. See the",
              f"separate [energy comparison]({energy_link}) and",
              f"[algorithm explanation]({algorithm_link}) for the accuracy cost.",
              "Timings are specific to this machine and workload. Plot error bars show",
              "ranges over repetitions, not confidence intervals.", "",
              "## Reproduce", "", "From the repository root:", "", "```bash",
              "python3 validation/benchmark_speed.py --repeats 3", "```", "",
              "Use `--eigen-include /path/to/eigen`, `--omp-threads 1 2 4`, or",
              "`--output-dir /path/to/report` as needed. Raw jobs are temporary.",
              "Full data: [timing_summary.csv](timing_summary.csv),",
              "[timing_samples.csv](timing_samples.csv), and [results.json](results.json).", ""]
    (output/"README.md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eigen-include", default=os.environ.get("EIGEN_INC", "/usr/include/eigen3"))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--omp-threads", nargs="+", type=int, default=[1, 2, 4])
    parser.add_argument("--output-dir", type=Path, default=ROOT/"validation/results/performance")
    args = parser.parse_args()
    if args.repeats < 3: parser.error("Use at least three timed repetitions")
    if not args.omp_threads or any(n <= 0 for n in args.omp_threads): parser.error("Thread counts must be positive")
    if not (Path(args.eigen_include)/"Eigen/Core").is_file(): parser.error("Eigen/Core not found")
    args.omp_threads = sorted(set(args.omp_threads))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    generator = load_generator()
    configurations = [(VARIANTS[0], 1), (VARIANTS[1], 1)] + [(VARIANTS[2], n) for n in args.omp_threads]
    rng = random.Random(20261021)
    samples, summaries, cases = [], [], []
    with tempfile.TemporaryDirectory(prefix="md-speed-") as temporary:
        work = Path(temporary)
        binaries = build(args.eigen_include, work)
        for grid, fraction, steps in CASES:
            count = grid**3
            case = f"N{count}_light{fraction:g}"
            job = work/case
            job.mkdir()
            cfg = dict(rc=2.5, sigma=1, epsilon=1, mass1=1, mass2=100, NPartRatio=fraction,
                       dt=0.01, T=2, totsteps=steps, nx=grid, ny=grid, nz=grid, rho=0.5)
            write_config(job, cfg)
            with contextlib.redirect_stdout(io.StringIO()):
                generator.generate_md_initial_conf(grid, grid, grid, .5, 2, 1, 1, 2.5, 1, 100,
                                                   fraction, filename=job/"initial_conf.csv", seed=42)
            case_info = {"case": case, "particles": count, "light_fraction": fraction, "macro_steps": steps,
                         "duration": steps*.01, "initial_csv_sha256": hashlib.sha256((job/"initial_conf.csv").read_bytes()).hexdigest(),
                         **pair_work(count, fraction)}
            cases.append(case_info)
            print(f"Case {case}: {steps} macro steps", flush=True)
            for mode in ("kernel", "application"):
                for variant, threads in configurations:
                    measure(mode, variant, threads, job, binaries)  # Discarded warmup.
                for repetition in range(args.repeats):
                    order = configurations.copy()
                    rng.shuffle(order)
                    for variant, threads in order:
                        seconds = measure(mode, variant, threads, job, binaries)
                        sample = {"case": case, "mode": mode, "variant": variant, "threads": threads,
                                  "repetition": repetition+1, "seconds": seconds}
                        samples.append(sample)
                for variant, threads in configurations:
                    times = [s["seconds"] for s in samples if s["case"]==case and s["mode"]==mode
                             and s["variant"]==variant and s["threads"]==threads]
                    summaries.append({"case": case, "particles": count, "light_fraction": fraction,
                                      "macro_steps": steps, "duration": steps*.01, "mode": mode,
                                      "variant": variant, "threads": threads,
                                      "median_seconds": statistics.median(times), "min_seconds": min(times),
                                      "max_seconds": max(times)})
                baseline = next(s["median_seconds"] for s in summaries if s["case"]==case and s["mode"]==mode and s["variant"]==VARIANTS[0])
                for summary in summaries:
                    if summary["case"]==case and summary["mode"]==mode:
                        summary["speedup_vs_uniform"] = baseline/summary["median_seconds"]
                        print(f"  {mode}: {summary['variant']} {summary['threads']}T: {summary['median_seconds']:.4f}s, {summary['speedup_vs_uniform']:.3f}x", flush=True)
    cpu_model = "unknown"
    for line in Path("/proc/cpuinfo").read_text().splitlines():
        if line.startswith("model name"):
            cpu_model = line.split(":", 1)[1].strip()
            break
    metadata = {"cpu_model": cpu_model, "logical_cpus": os.cpu_count(),
                "allowed_cpus": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
                "platform": platform.platform(), "compiler": run([*shlex.split(os.environ.get("CXX", "g++")), "--version"]).stdout.splitlines()[0],
                "build_mode": "releaseO3", "repeats": args.repeats, "warmups_per_case_mode_configuration": 1,
                "omp_settings": {k: environment(1)[k] for k in ("OMP_DYNAMIC", "OMP_PROC_BIND", "OMP_PLACES", "OMP_WAIT_POLICY")},
                "dt": .01, "density": .5, "masses": [1, 100], "ratio": 10, "seed": 42,
                "fairness": "Same initial CSV, macro step count and physical duration within each case; not matched-error timing.",
                "scope": "One machine, warm filesystem caches, sequential runs; no claim of exclusive CPU access or universal speedup."}
    result = {"metadata": metadata, "cases": cases, "summary": summaries, "samples": samples}
    (args.output_dir/"results.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    write_csv(args.output_dir/"timing_summary.csv", summaries)
    write_csv(args.output_dir/"timing_samples.csv", samples)
    plot_speed(summaries, args.output_dir, configurations)
    report(args.output_dir, metadata, summaries, cases)
    print(f"Saved benchmark to {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
