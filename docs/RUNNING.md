# Building, running and configuring simulations

These instructions apply to all three programs. For the update schedule and
variant-specific implementation, use the guides for
[Normal](../MDNormal/docs/README.md),
[serial MTS](../MDMultipleTimeSteps/docs/README.md), and
[OpenMP MTS](../MDMultipleTimeStepsOpenMP/docs/README.md).
Commands inside a variant directory are identified below.

## Build

Requirements:

- GNU Make and a C++17 compiler. Local verification used GCC 7.5.
- Eigen headers. Local verification used Eigen 3.4.90 headers (the development snapshot available locally).
- GCC's OpenMP runtime for `MDMultipleTimeStepsOpenMP`.
- `libstdc++fs`, because the I/O wrappers use `std::experimental::filesystem`.
- Python 3.8+ for initial-condition generation and the regression suite.

On Debian/Ubuntu, the build dependencies are available as:

```bash
sudo apt install build-essential libeigen3-dev python3
```

Build from a variant directory:

```bash
cd MDNormal
make -j2 BUILD_MODE=releaseO3
```

Eigen defaults to `/usr/include/eigen3`. For another installation, pass the
directory that **contains** `Eigen/`:

```bash
make -j2 BUILD_MODE=releaseO3 EIGEN_INC=/path/to/eigen
```

| Mode | Flags | Executable |
| :--- | :--- | :---: |
| `debug` (default) | `-g -O0` | `bin/MolecularDynamics_debug` |
| `releaseO2` | `-g -O2` | `bin/MolecularDynamics_O2` |
| `releaseO3` | `-g -O3` | `bin/MolecularDynamics_O3` |

Objects are separated by build mode, and header dependencies are tracked.
Run `make clean` before changing compiler, compiler flags, or Eigen installation.
It removes generated builds and all three executable modes.

## Run

Each variant includes a ready-to-run `jobs/sim1/` with a configuration and
particle CSV. From the same variant directory used to build:

```bash
./bin/MolecularDynamics_O3 sim1
```

For OpenMP:

```bash
cd MDMultipleTimeStepsOpenMP  # from the repository root
make -j2 BUILD_MODE=releaseO3
OMP_NUM_THREADS=4 ./bin/MolecularDynamics_O3 sim1
```

A job name resolves to `jobs/<name>/` beside the executable's `bin/` directory.
You may instead pass an existing job directory, including an absolute path:

```bash
./bin/MolecularDynamics_O3 /path/to/my-job
```

**Running a job overwrites its energy log and snapshots with matching step
numbers.** Copy the job to a new directory to preserve a previous run. Snapshots
from other step numbers remain in the directory; use a fresh job directory
when comparing runs.

### Generate initial conditions

After changing grid size, density, masses, fraction, or initial temperature,
regenerate the particle CSV to match the configuration:

```bash
python3 jobs/sim1/generateIC.py --seed 42
```

The script defaults to the config beside itself. It also supports explicit
paths, useful for a new job:

```bash
python3 jobs/sim1/generateIC.py --config jobs/my-job/config.txt --seed 42
```

This writes `initial_conf.csv` beside the selected config; `--output` overrides
the destination. The default seed is 0. Positions form a centered simple-cubic
lattice. Velocities are Gaussian with standard deviation `sqrt(T/m)`, using
`k_B = 1`, followed by subtraction of center-of-mass velocity. The sampled
velocities are not rescaled to exactly match `T`, and `T` is not maintained
during simulation.

The input CSV must have the header and column order:

```csv
id,m,x,y,z,vx,vy,vz
```

IDs must be unique, fields finite, and the particle count equal to
`nx * ny * nz`. Masses must match the two configured values. Periodic overlaps
are rejected. `NPartRatio` determines the generated mixture; when supplying a
custom CSV, the actual mixture comes from the CSV.

## Configuration

All keys below are required. Use `key value` or `key = value`. Blank lines and
`#` comments are supported. Unknown, duplicate, missing, and invalid values
produce an error and a nonzero exit status. The old `delimiter = 0` line is
accepted for compatibility; particle CSV always uses commas.

| Key | Meaning and constraints |
| --- | --- |
| `nx`, `ny`, `nz` | Positive integer lattice counts; their product is the particle count |
| `rho` | Positive number density, N/volume; lattice spacing is `rho^(-1/3)` |
| `sigma` | Positive LJ zero-crossing distance |
| `epsilon` | Positive LJ well depth |
| `rc` | Positive cutoff, at most half the shortest box dimension |
| `mass1`, `mass2` | Positive masses; their order does not matter |
| `NPartRatio` | Fraction of the **lighter** species generated, in [0, 1] |
| `T` | Nonnegative initial velocity sampling temperature, with `k_B = 1` |
| `dt` | Positive macro time step; also the interval between energy rows |
| `totsteps` | Nonnegative integer number of macro steps |

The box dimensions are `(nx, ny, nz) * rho^(-1/3)`. All lengths, masses,
energies, and times must use a consistent unit system; the examples use
`sigma = epsilon = m_light = 1`.

The pair potential inside the cutoff is
`U(r) = 4*epsilon*((sigma/r)^12 - (sigma/r)^6) - U(rc)`, and zero outside.
This constant shift makes the potential continuous at `rc`; the force is
truncated without a force shift and is discontinuous there. There is no
long-range tail correction. Distances use the minimum image convention.

## Output

All output is written into the selected job directory. `OutputTimes/` is
created automatically.

- `energy.csv`: columns `t,K,U,Etot`. `t` is the **integer macro step index**;
  physical time is `t*dt`. Includes the initial and final states, giving
  `totsteps + 1` rows. `U` counts each pair once.
- `OutputTimes/OutputTimes<step>.txt`: comma-separated snapshots at step 0,
  every 100 steps, and the final step, with columns
  `id,m,x,y,z,vx,vy,vz,ax,ay,az,fx,fy,fz,fxAux,fyAux,fzAux,U`.
  Accelerations, forces, and per-particle potential energy are recomputed on
  copies at the snapshot positions. Each particle receives half of each pair
  energy, so its `U` column sums to the total potential. In MTS, `fAux`
  contains the integrator's saved same-group force buffer; it is a diagnostic, not an
  additional physical force to add to the snapshot's `f`.

Numeric output uses round-trip double precision. Non-finite energy or particle
state stops the run with a diagnostic; reduce `dt` and inspect the inputs.
The final integration step is checked and logged as well.
