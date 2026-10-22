# OpenMP Tuckerman disparate-mass MTS

This program parallelizes the same numerical update as
[serial MTS](../../MDMultipleTimeSteps/docs/README.md). It uses the same
light/heavy/light schedule, fixed ratio, restricted Verlet steps, force
caches and energy bookkeeping. Threading does not define a different
integrator. The [algorithm guide](../../docs/ALGORITHM.md) explains the
method and its comparison with uniform small-step Verlet.

## Build and run

From the repository root:

```bash
make -C MDMultipleTimeStepsOpenMP -j2 BUILD_MODE=releaseO3
OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE ./MDMultipleTimeStepsOpenMP/bin/MolecularDynamics_O3 sim1
```

The example is `MDMultipleTimeStepsOpenMP/jobs/sim1/`. Set `OMP_NUM_THREADS`
to the desired team size. GCC's OpenMP runtime is required. Eigen defaults
to `/usr/include/eigen3`; pass `EIGEN_INC=/path/to/eigen` to `make` if needed.

Running a job overwrites its energy log and matching snapshots. Copy the job
to preserve previous output. See the [shared run guide](../../docs/RUNNING.md)
for dependencies, configuration, generation and output schemas. As in serial
MTS, both mass groups must be distinct and nonempty.

## Parallel implementation

The core is
[`MDParticleEvolver.cpp`](../src/MDParticleEvolver/MDParticleEvolver.cpp).
The serial guide's [force-buffer table](../../MDMultipleTimeSteps/docs/README.md#force-buffers)
also applies here.

| Operation | Parallel mechanism |
| --- | --- |
| Same-group and cross-group pair loops | Work shared among threads |
| Shared particle force components | Atomic additions |
| Pair potential and kinetic energy | Reductions |
| Acceleration refresh | After force accumulation completes |
| Kicks, drifts and buffer operations | Parallel loops over distinct particles |

The ends of the parallel loops provide the barriers needed before the next
stage consumes the forces. Atomic force additions and energy reductions do
not impose a unique floating-point summation order.

## Numerical agreement and speed

One-thread OpenMP is the control for the serial update. Multiple threads
can introduce roundoff-scale force differences that grow along a long
trajectory. Two individual maximum energy errors are not a systematic
parallel accuracy ratio. The optional
[parallel report](../../validation/results/parallel_accuracy/README.md)
checks identical states, energy sums and sensitivity to tiny input changes.

Threading also adds overhead from team setup, synchronization and atomic
updates. A small job can benefit little or run slower with more threads.
Use the [speed report](../../validation/results/performance/README.md) for
measured medians, ranges, particle counts and team settings. It distinguishes
integration-only time from complete runtime, which includes output.

The [energy report](../../validation/results/energy_order/README.md) compares
the implementations on identical inputs. The
[algorithm guide](../../docs/ALGORITHM.md) explains why MTS's coarser heavy
steps and frozen-group splitting can produce a larger error amplitude than
the uniform baseline, even when both are second order.

## Validation

Run `python3 validation/run_tests.py --build-mode releaseO3` from the
repository root. The suite checks serial/OpenMP agreement and the uncached
Tuckerman reference. All scripts and reports remain under the single
[validation folder](../../validation/README.md).
