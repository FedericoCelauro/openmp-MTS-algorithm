# Tests and validation

There are three main commands, all run from the repository root:

| Goal | Command | Output |
| --- | --- | --- |
| Check implementation correctness | `python3 validation/run_tests.py --build-mode releaseO3` | Pass/fail in the terminal |
| Compare energy errors and integration order | `python3 validation/check_energy_order.py` | [Energy report](results/energy_order/README.md) |
| Compare execution speed | `python3 validation/benchmark_speed.py --repeats 3` | [Speed report](results/performance/README.md) |

Use `--eigen-include /path/to/eigen` or set `EIGEN_INC` if Eigen is not
installed under `/usr/include/eigen3`. The C++ requirements are in the
[shared run guide](../docs/RUNNING.md#build).

The correctness suite needs only Python's standard library. The energy
and speed scripts also need:

```bash
python3 -m pip install -r requirements-analysis.txt
```

All runs use temporary jobs. The analysis scripts save CSVs, plots and a
report under `validation/results/`; `--output-dir` selects another location.
They do not depend on precomputed results. Run timing benchmarks on their
own, without another validation run competing for the CPU.

## What does each command establish?

`run_tests.py` has 15 test groups covering pair forces, periodic boundaries,
input/output, energy reconstruction, serial/OpenMP agreement, reversibility,
smooth-force convergence, and agreement with the Tuckerman update. Its
independent reference recalculates all forces before every kick, so it also
checks the production force caches. Both MTS variants are checked on four
interacting particles for ratios 2 and 10, plus the odd-ratio extensions at
1 and 3. This is a correctness check for the exercised cases, not a universal
stability guarantee.

`check_energy_order.py` starts every run from the same particle CSV and
configuration. It compares all three implementations for physical duration
50 at macro steps 0.04, 0.05, 0.08 and 0.1. It includes one-thread OpenMP
as a serial control. The report shows energy error relative to the initial
value and fits the standard deviation over the second half of each run.
Options `--dt`, `--duration` and `--threads` change the experiment.

`benchmark_speed.py` uses identical inputs within each case, an optimized
build, a discarded warmup and at least three timed repetitions. It reports
medians and ranges for integration alone and the complete application.
Options `--repeats` and `--omp-threads` select the repetitions and thread
counts. The comparison uses equal macro steps and durations, rather than
matched energy error.

For the reason MTS has a larger energy error at the same macro step, see
[the algorithm guide](../docs/ALGORITHM.md#why-can-mts-have-a-larger-energy-error).

## Folder layout

```text
validation/
  README.md                   start here
  run_tests.py                implementation checks
  check_energy_order.py       accuracy comparison and energy plots
  benchmark_speed.py          speed comparison
  support/                    reference code and compiled helper drivers
  results/                    generated reports, figures and numerical data
```

The files in `support/` are used by the main commands. They are not separate
tests that must each be run. The optional diagnostic
`python3 validation/support/check_parallel_accuracy.py` investigates
serial/parallel trajectory sensitivity; its recorded
[report](results/parallel_accuracy/README.md) is separate from the speed test.

## Other build modes and sanitizers

```bash
python3 validation/run_tests.py --build-mode debug
python3 validation/run_tests.py --build-mode releaseO2
python3 validation/run_tests.py --sanitize
```

The sanitizer command builds temporary AddressSanitizer/UndefinedBehaviorSanitizer
binaries. If LeakSanitizer cannot run under a debugger or tracer, use
`ASAN_OPTIONS=detect_leaks=0` for that local run; leak checking is then disabled.
Continuous integration runs all three build modes and a sanitizer job.
