# Tuckerman Disparate-Mass Molecular Dynamics

C++17 reproduction of the disparate-mass `G_xyx` algorithm in Section IV.A
of Tuckerman, Berne and Martyna (1992),
[*Reversible multiple time scale molecular dynamics*](https://doi.org/10.1063/1.463137),
*Journal of Chemical Physics* **97**, 1990–2001.
The implementation uses the velocity-Verlet factors in equations (4.2)–(4.6).

Three programs compare this method with uniform small-step Verlet for
Lennard-Jones particles with periodic boundaries. The numerical experiments
use this project's inputs rather than reproduce every figure in the paper.

## Programs

| Program | Integration and force calculation | Guide |
| --- | --- | --- |
| `MDNormal` | Every particle takes small Verlet steps; serial forces | [Normal docs](MDNormal/docs/README.md) |
| `MDMultipleTimeSteps` | Light/heavy/light MTS; serial forces | [Serial MTS docs](MDMultipleTimeSteps/docs/README.md) |
| `MDMultipleTimeStepsOpenMP` | Same MTS update; parallel forces | [OpenMP MTS docs](MDMultipleTimeStepsOpenMP/docs/README.md) |

The small step is `dt/Ratio`, where `Ratio` is fixed from the mass ratio.
MTS uses the full `dt` for heavy motion and holds the inactive group's
positions and velocities fixed. See the [algorithm guide](docs/ALGORITHM.md)
for the schedule, odd-ratio extension and energy-error explanation.

## Quick start

Requires GNU Make, a C++17 compiler and Eigen headers. The parallel program
also requires OpenMP. From the repository root:

```bash
make -C MDNormal -j2 BUILD_MODE=releaseO3
./MDNormal/bin/MolecularDynamics_O3 sim1
```

Each program includes a runnable `jobs/sim1/`. Running it overwrites its
energy log and matching snapshots. Build/run commands for the other programs
are in their guides above. For installation, other build modes, configuration
keys, initial-condition generation and output columns, use
[Building and running](docs/RUNNING.md).

## Compare accuracy and speed

All checks, analysis scripts and reports are under
[validation/](validation/README.md). From the repository root:

```bash
python3 validation/run_tests.py --build-mode releaseO3
python3 validation/check_energy_order.py
python3 validation/benchmark_speed.py --repeats 3
```

The first command needs only Python's standard library; the last two also
need `python3 -m pip install -r requirements-analysis.txt`. Supply
`--eigen-include /path/to/eigen` or `EIGEN_INC` for a custom Eigen installation.

The [energy report](validation/results/energy_order/README.md) compares
common initial states, macro steps and durations. The
[speed report](validation/results/performance/README.md) gives integration
and complete-application timings with repetitions and thread settings.

In the controlled energy runs, MTS has larger error amplitudes at the same
macro step: the normal baseline also gives heavy particles the smaller step.
This is not a universal accuracy ranking. Likewise, fixed-step speedups are
not matched-accuracy gains. The algorithm guide explains these distinctions.

## Documentation

- [Building and running](docs/RUNNING.md): shared dependencies, inputs and outputs.
- [Algorithm comparison](docs/ALGORITHM.md): update schedules and energy errors.
- [Source overview](docs/IMPLEMENTATION.md): shared classes and per-program implementation guides.
- [Validation](validation/README.md): the three commands, results and optional diagnostics.

Copyright (c) 2026 Federico Celauro. Released under the [MIT License](LICENSE).
