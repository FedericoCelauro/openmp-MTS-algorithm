# Source overview

The programs use separate source trees so the update schedules and force
loops can be compared directly. Each has `src/`, `src/include/`, a Makefile
and runnable job inputs. The detailed implementation lives in each
program's own documentation:

| Program | Implementation guide |
| --- | --- |
| Uniform small-step Verlet | [MDNormal](../MDNormal/docs/README.md) |
| Serial Tuckerman MTS | [MDMultipleTimeSteps](../MDMultipleTimeSteps/docs/README.md) |
| OpenMP Tuckerman MTS | [MDMultipleTimeStepsOpenMP](../MDMultipleTimeStepsOpenMP/docs/README.md) |

## Shared classes

```text
PointInSpace       particle ID and position
  LagrangianPart   velocity and position updates
    MassivePart   mass, acceleration and velocity kicks
      LJPart      interaction parameters and force storage
        MDPart    periodic geometry and cutoff
```

`SimParams` validates configuration and derives box dimensions and the
fixed step ratio. `MDParticleEvolver` owns the particles, force calculations
and integration. `MDTimeIntegrator` schedules macro steps, energy logging
and snapshots. `MolecularDynamics` resolves jobs and reports startup errors.
The ASCII wrappers manage the input/output streams.

The same pair potential, input format and output conventions are used by all
three programs; see [RUNNING.md](RUNNING.md). Shared force-law and I/O fixes
should be applied consistently. The single [validation suite](../validation/README.md)
checks all three, including the two MTS implementations against an uncached
reference to catch mistakes shared by their optimized buffers.
