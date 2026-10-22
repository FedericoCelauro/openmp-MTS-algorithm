# Uniform small-step velocity Verlet

`MDNormal` is the serial reference implementation. Every particle takes
`Ratio` velocity-Verlet steps of duration `dt/Ratio` during each macro step.
For masses 1 and 100, `Ratio = 10`: both light and heavy particles use the
small step. This is the baseline used by the energy and speed comparisons.
It is not a single full-system Verlet step of duration `dt`.

## Build and run

From the repository root:

```bash
make -C MDNormal -j2 BUILD_MODE=releaseO3
./MDNormal/bin/MolecularDynamics_O3 sim1
```

The example is `MDNormal/jobs/sim1/`. A job name is resolved relative to
the program's project directory; an existing job-directory path also works.
Running a job overwrites its energy log and matching snapshots. Use a copy
of the job directory to preserve earlier output.

Eigen defaults to `/usr/include/eigen3`; pass `EIGEN_INC=/path/to/eigen`
to `make` for another installation. See the
[shared run guide](../../docs/RUNNING.md) for dependencies, build modes,
initial-condition generation, configuration keys and output columns.

## Update and force calculation

The entry points are in
[`MDParticleEvolver.cpp`](../src/MDParticleEvolver/MDParticleEvolver.cpp).

| Function | Role |
| --- | --- |
| `evolvePriv` | Repeat the small step `Ratio` times |
| `VelocityVerlet` | Half-kick every particle, drift every particle, recompute forces, half-kick every particle |
| `calcForcesPairwise` | Visit unordered pairs over the full particle collection |
| `writeStatePriv` | Write physical particle diagnostics |

Each pair force contributes equally and oppositely to the two particles.
All momenta are kicked in the same full-system stages. The next step reuses
the preceding endpoint forces. The pair potential, periodic boundaries and
cutoff are the same as in both MTS variants.

`MDNormal` supports a single species as well as two mass groups. The time
ratio is fixed from the configured masses, with a minimum of one.

## Accuracy and computational cost

Full-system Verlet is reversible and canonical symplectic for smooth
Hamiltonian forces. Its physical energy is approximately conserved rather
than exactly constant at finite step. The force discontinuity at the cutoff
limits the usual smooth-force energy and convergence arguments.

At a common macro `dt`, this baseline resolves heavy motion more frequently
than MTS. In the controlled comparison it has smaller energy fluctuations,
at the cost of recalculating all candidate pairs at every small step.
Different methods need not have equal error coefficients even when they
are both second order. See the
[algorithm comparison](../../docs/ALGORITHM.md) for a schedule and energy plot.

The [energy report](../../validation/results/energy_order/README.md) and
[speed report](../../validation/results/performance/README.md) use common
inputs within each comparison. Speedups at equal macro step are not
matched-accuracy speedups.

## Validation

From the repository root:

```bash
python3 validation/run_tests.py --build-mode releaseO3
```

This checks all three programs. For the energy and speed commands, see
the single [validation guide](../../validation/README.md).
