"""Generate a simple-cubic LJ configuration with Maxwell-Boltzmann velocities.

Units use k_B = 1. The center-of-mass velocity is removed; temperature is
sampled, not rescaled to an exact value. NPartRatio refers to the lighter mass.
"""

import argparse
import csv
import math
from pathlib import Path
import random


KEYS = {"nx", "ny", "nz", "rho", "T", "sigma", "epsilon", "rc", "mass1",
        "mass2", "NPartRatio", "dt", "totsteps"}
INTEGER_KEYS = {"nx", "ny", "nz", "totsteps"}


def load_config(filename):
    """Read the same required key/value format as the C++ configuration parser."""
    config = {}
    for original in Path(filename).read_text().splitlines():
        line = original.split("#", 1)[0].replace("=", " ").split()
        if not line:
            continue
        if len(line) != 2:
            raise ValueError(f"Expected a key and value: {original}")
        key, value = line
        if key == "delimiter" and value == "0":  # Legacy unused setting.
            continue
        if key not in KEYS or key in config:
            raise ValueError(f"Unknown or duplicate config key: {key}")
        config[key] = int(value) if key in INTEGER_KEYS else float(value)
        if not math.isfinite(config[key]):
            raise ValueError(f"Non-finite value for {key}")
    missing = KEYS - config.keys()
    if missing:
        raise ValueError(f"Missing config keys: {', '.join(sorted(missing))}")
    if config["dt"] <= 0 or config["totsteps"] < 0:
        raise ValueError("dt must be positive and totsteps nonnegative")
    return config


def generate_md_initial_conf(nx, ny, nz, rho, T, sigma, epsilon, rc,
                             mass1, mass2, n_ratio, filename="initial_conf.csv",
                             seed=0):
    """Write eight CSV columns; the first int(N*n_ratio) particles are light.

    A local random generator makes the default seed reproducible and leaves
    Python's global random state unchanged. Use a different seed for new samples.
    """
    if any(not isinstance(n, int) or n <= 0 for n in (nx, ny, nz)):
        raise ValueError("Grid dimensions must be positive integers")
    if any(not math.isfinite(x) or x <= 0 for x in
           (rho, sigma, epsilon, rc, mass1, mass2)):
        raise ValueError("Density, LJ parameters, cutoff and masses must be positive and finite")
    if not math.isfinite(T) or T < 0 or not math.isfinite(n_ratio) or not 0 <= n_ratio <= 1:
        raise ValueError("Require finite T >= 0 and NPartRatio in [0, 1]")

    count = nx * ny * nz
    spacing = rho ** (-1.0 / 3.0)  # rho is an absolute number density.
    box = (nx * spacing, ny * spacing, nz * spacing)
    if rc > min(box) / 2:
        raise ValueError("rc must be at most half the shortest box dimension")
    light, heavy = sorted((mass1, mass2))
    n_light = int(count * n_ratio)
    rng = random.Random(seed)
    rows = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                mass = light if len(rows) < n_light else heavy
                pos = [(i + 0.5) * spacing - length / 2
                       for i, length in zip((ix, iy, iz), box)]
                vel = [math.sqrt(T / mass) * rng.gauss(0, 1) for _ in range(3)]
                rows.append([len(rows), mass, *pos, *vel])
    total_mass = sum(row[1] for row in rows)
    v_cm = [sum(row[1] * row[5 + k] for row in rows) / total_mass for k in range(3)]
    for row in rows:
        for k in range(3):
            row[5 + k] -= v_cm[k]
    with open(filename, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["id", "m", "x", "y", "z", "vx", "vy", "vz"])
        writer.writerows(rows)
    print(f"Generated {count} particles in {filename}; seed={seed}, box={box}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path(__file__).with_name("config.txt"))
    parser.add_argument("--output", type=Path, help="Defaults to initial_conf.csv beside the config")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    cfg = load_config(args.config)
    generate_md_initial_conf(
        **{key: cfg[key] for key in ("nx", "ny", "nz", "rho", "T", "sigma", "epsilon", "rc", "mass1", "mass2")},
        n_ratio=cfg["NPartRatio"], filename=args.output or args.config.with_name("initial_conf.csv"),
        seed=args.seed,
    )
