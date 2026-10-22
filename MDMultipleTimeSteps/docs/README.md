# Serial Tuckerman disparate-mass MTS

This program implements the velocity-Verlet realization of Tuckerman,
Berne and Martyna's `G_xyx` algorithm, Section IV.A, equations (4.2)–(4.6).
The [algorithm guide](../../docs/ALGORITHM.md) gives the reference, update
schedule and energy-error comparison with uniform small-step Verlet.

## Build and run

From the repository root:

```bash
make -C MDMultipleTimeSteps -j2 BUILD_MODE=releaseO3
./MDMultipleTimeSteps/bin/MolecularDynamics_O3 sim1
```

The example is `MDMultipleTimeSteps/jobs/sim1/`. Running a job overwrites
its energy log and matching snapshots; copy the job directory to preserve
earlier output. Eigen defaults to `/usr/include/eigen3`; use
`EIGEN_INC=/path/to/eigen` with `make` for another installation.
The [shared run guide](../../docs/RUNNING.md) documents dependencies,
configuration, initial-condition generation and output columns.

MTS requires two distinct, nonempty mass groups. They are selected by actual
mass, independently of the order of `mass1` and `mass2`. For a single species,
use [MDNormal](../../MDNormal/docs/README.md).

## Stepper and source map

The core is
[`MDParticleEvolver.cpp`](../src/MDParticleEvolver/MDParticleEvolver.cpp).

| Function | Role |
| --- | --- |
| `evolvePriv` | Light half-sequence, heavy full step, light half-sequence |
| `VelocityVerlet` | Kick and drift only the active group; recompute its forces |
| `calcForcesPairwiseSame` | Same-group forces and their saved auxiliary buffer |
| `calcForcesPairwiseDiff` | Cross forces on both groups |
| `combineForces` | Restore the complete force when switching active group |
| `writeStatePriv` | Recompute physical snapshot diagnostics on copies |

The fixed ratio is `max(1, floor(sqrt(m_heavy/m_light)))`. For an even ratio,
each light half-sequence has `Ratio/2` steps of `dt/Ratio`. Odd ratios use a
mirrored half-small-step on each side so that every light stage still has
duration `dt/2`. That odd-ratio handling is this project's extension.

Both positions and velocities of the inactive group stay fixed. Cross forces
are nevertheless calculated on it for later use. This is the selected
coordinate-group splitting; adding inactive momentum kicks would change
the algorithm.

## Force buffers

`LJPart::setForce` stores `f` and acceleration `f/m`; it does not change
velocity. `MassivePart::updateVel(s)` performs the velocity kick.

- `fAux` saves the latest same-group force.
- `f` is the working force. It contains either the complete active force
  or the inactive group's cross force, depending on the stage.

The same-group routine clears the active working force and saves its result
in `fAux`. The cross routine clears the inactive working force, adds cross
forces to the active one, and stores the opposite cross forces on the
inactive group. This refreshes a force sample; it does not accumulate an
impulse over the light substeps.

Write `LL`, `HH`, `LH`, and `HL` for force contributions, where `LH` acts
on light particles and `HL` acts on heavy particles:

| Stage boundary | Light `f` | Heavy `f` | Saved same-group buffers |
| --- | --- | --- | --- |
| After a light substep | `LL + LH` | `HL` | Light `LL`, heavy `HH` |
| After `combineForces(heavy)` | `LL + LH` | `HH + HL` | Unchanged |
| After the heavy substep | `LH` | `HH + HL` | Heavy `HH` refreshed |
| After `combineForces(light)` | `LL + LH` | `HH + HL` | Unchanged |
| After the final light substep | `LL + LH` | `HL` | Light `LL` refreshed |

Initialization computes complete forces for both groups. The first light
substep establishes the first table row before the first group switch.
Heavy–heavy forces remain valid while only light positions move, and
light–light forces remain valid while only heavy positions move. Cross
forces are refreshed whenever either group moves. The initial kick of a
restricted Verlet reuses the complete active force at its starting state.

## Potential energy and snapshots

After a light substep, working `U` contains `V_LL + V_LH`. The heavy step
saves its recalculated `V_HH` in `UAux`. Heavy positions stay fixed through
the last light sequence, so that contribution remains valid. At the end of
the macro step, `U += UAux` gives the physical total potential, counting
each pair once. Initialization also counts all pair contributions once.

Energy is logged at macro-step boundaries. Intermediate buffers are working
data rather than complete physical diagnostics. Snapshots recompute full
forces, accelerations and half-pair energies on copies so they cannot alter
the next integration kick. Snapshot `f` is therefore a physical force;
snapshot `fAux` is the saved same-group buffer, not an extra force to add.

## Accuracy, speed and validation

The selected fixed-step method is reversible and second order for smooth
forces, but not generally canonical symplectic. Pair forces are opposite,
while their momentum kicks occur in different stages; exact full-system
momentum conservation does not follow from the pair-force symmetry alone.

MTS saves heavy–heavy force work and can have larger energy errors than
uniform small-step Verlet at the same macro step. The
[energy](../../validation/results/energy_order/README.md) and
[speed](../../validation/results/performance/README.md) reports quantify the
tested cases without claiming a universal accuracy ranking or speedup.

All validation stays in one folder. Run
`python3 validation/run_tests.py --build-mode releaseO3` from the repository
root, or see the [validation guide](../../validation/README.md) for the
energy and speed comparisons. The suite compares the optimized buffers with
an [uncached reference](../../validation/support/tuckerman_reference.py)
that recalculates every pair force before each kick.
