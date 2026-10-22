"""Uncached reference for the disparate-mass G_xyx update.

See Tuckerman, Berne and Martyna, J. Chem. Phys. 97, 1990 (1992),
Section IV.A, Eqs. (4.2)-(4.6), DOI: 10.1063/1.463137.

This deliberately recomputes every pair force before every kick and uses no
engine classes or force caches. It is a small-system numerical reference,
not a performance implementation. Odd ratios use this project's mirrored
half-small-step extension to the paper's even-ratio schedule.
"""

import math


def advance(initial, config):
    """Return positions and velocities after config['totsteps'] macro steps.

    initial contains rows [mass, x, y, z, vx, vy, vz]. The force law is the
    project's truncated Lennard-Jones force; the constant potential shift
    does not affect it. Mass groups and step ratio remain fixed throughout.
    """
    masses = [row[0] for row in initial]
    positions = [list(row[1:4]) for row in initial]
    velocities = [list(row[4:7]) for row in initial]
    light_mass, heavy_mass = min(masses), max(masses)
    if light_mass == heavy_mass:
        raise ValueError("The disparate-mass reference needs two mass groups")
    light = [i for i, mass in enumerate(masses) if mass == light_mass]
    heavy = [i for i, mass in enumerate(masses) if mass == heavy_mass]
    if len(light) + len(heavy) != len(masses):
        raise ValueError("Exactly two mass groups are supported")
    spacing = config["rho"] ** (-1.0 / 3.0)
    box = [config[key] * spacing for key in ("nx", "ny", "nz")]
    cutoff2 = config["rc"] ** 2
    sigma2 = config["sigma"] ** 2
    epsilon = config["epsilon"]

    def forces():
        result = [[0.0, 0.0, 0.0] for _ in positions]
        for i in range(len(positions)):
            for j in range(i):
                displacement = [positions[i][d] - positions[j][d] for d in range(3)]
                displacement = [value - length * round(value / length)
                                for value, length in zip(displacement, box)]
                distance2 = sum(value * value for value in displacement)
                if distance2 == 0:
                    raise ValueError("Overlapping reference particles")
                if distance2 <= cutoff2:
                    scaled6 = (sigma2 / distance2) ** 3
                    factor = 24 * epsilon * (2 * scaled6 * scaled6 - scaled6) / distance2
                    for d in range(3):
                        pair_force = factor * displacement[d]
                        result[i][d] += pair_force
                        result[j][d] -= pair_force
        return result

    def kick(active, duration):
        current = forces()
        for i in active:
            for d in range(3):
                velocities[i][d] += duration * current[i][d] / masses[i]

    def verlet(active, duration):
        kick(active, duration / 2)
        for i in active:
            for d in range(3):
                positions[i][d] += duration * velocities[i][d]
                positions[i][d] -= box[d] * round(positions[i][d] / box[d])
        kick(active, duration / 2)

    ratio = max(1, int(math.sqrt(heavy_mass / light_mass)))
    macro = config["dt"]
    first_half = [macro / ratio] * (ratio // 2)
    if ratio % 2:
        first_half.append(macro / (2 * ratio))
    for _ in range(config["totsteps"]):
        for duration in first_half:
            verlet(light, duration)
        verlet(heavy, macro)
        for duration in reversed(first_half):
            verlet(light, duration)
    return positions, velocities
