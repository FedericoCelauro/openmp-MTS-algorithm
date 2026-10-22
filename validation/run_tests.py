#!/usr/bin/env python3
"""Build all variants and run physics, I/O, input and integration regressions.

Only Python's standard library is required. Jobs and logs live in a temporary
directory, so this suite never replaces checked-in simulation data.
"""

import argparse
import csv
import importlib.util
import math
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

from support.tuckerman_reference import advance as tuckerman_reference

ROOT = Path(__file__).resolve().parents[1]
VARIANTS = ("MDNormal", "MDMultipleTimeSteps", "MDMultipleTimeStepsOpenMP")
CONFIG = {"rc": 2.5, "sigma": 1.0, "epsilon": 1.0, "mass1": 1.0,
          "mass2": 100.0, "NPartRatio": 0.25, "dt": 0.01, "T": 0.2,
          "totsteps": 20, "nx": 4, "ny": 4, "nz": 4, "rho": 0.5}
SPEC = importlib.util.spec_from_file_location("generate_ic", ROOT / VARIANTS[0] / "jobs/sim1/generateIC.py")
GENERATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GENERATOR)


def command(args, **kwargs):
    result = subprocess.run(list(map(str, args)), text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, timeout=180, **kwargs)
    if result.returncode:
        raise AssertionError(f"Command failed: {args}\n{result.stdout}\n{result.stderr}")
    return result


def write_config(path, cfg):
    path.write_text("# Regression fixture\n" + "\n".join(f"{k} = {v}" for k, v in cfg.items()) + "\n")


def fixture(path, **overrides):
    path.mkdir(parents=True, exist_ok=True)
    cfg = dict(CONFIG, **overrides)
    write_config(path / "config.txt", cfg)
    # Suppress routine generation output while preserving errors.
    command([os.sys.executable, ROOT / VARIANTS[0] / "jobs/sim1/generateIC.py",
             "--config", path / "config.txt", "--seed", "42"])
    return cfg


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def run_sim(variant, job, threads=1, succeeds=True):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE")
    result = subprocess.run([str(BINARIES[variant]), str(job)], cwd=job.parent,
                            env=env, capture_output=True, text=True, timeout=90)
    if succeeds:
        if result.returncode:
            raise AssertionError(result.stdout + result.stderr)
    else:
        if result.returncode <= 0:
            raise AssertionError("Invalid input must fail cleanly, without aborting\n" + result.stderr)
        if "Simulation failed:" not in result.stderr and "Simulation crashed:" not in result.stderr:
            raise AssertionError("Missing useful diagnostic: " + result.stderr)
    return result


def shifted_energy(snapshot, box, rc):
    """Independent scalar LJ calculation on serialized positions."""
    total = 0.0
    for i, a in enumerate(snapshot):
        for b in snapshot[i+1:]:
            r2 = 0.0
            for axis, length in zip("xyz", box):
                r = float(a[axis]) - float(b[axis])
                r -= length * round(r / length)
                r2 += r*r
            if r2 <= rc*rc:
                q6 = r2 ** -3
                total += 4 * (q6*q6 - q6 - rc**-12 + rc**-6)
    return total


class RegressionChecks(unittest.TestCase):
    def job(self, name, **overrides):
        path = TEMP / self._testMethodName / name
        cfg = fixture(path, **overrides)
        return path, cfg

    def test_physics_and_file_utilities(self):
        for variant in VARIANTS:
            with self.subTest(variant=variant):
                out = TEMP / "unit" / variant
                out.mkdir(parents=True)
                command([UNIT_BINARIES[variant], out])

    def test_sample_inputs_match_configuration(self):
        for variant in VARIANTS:
            with self.subTest(variant=variant):
                cfg = GENERATOR.load_config(ROOT / variant / "jobs/sim1/config.txt")
                data = rows(ROOT / variant / "jobs/sim1/initial_conf.csv")
                self.assertEqual(len(data), cfg["nx"] * cfg["ny"] * cfg["nz"])
                self.assertEqual(sum(float(p["m"]) == min(cfg["mass1"], cfg["mass2"]) for p in data),
                                 int(len(data) * cfg["NPartRatio"]))

    def test_generator_reproducible_and_zero_momentum(self):
        job, _ = self.job("generator")
        before = (job / "initial_conf.csv").read_bytes()
        command([os.sys.executable, ROOT / VARIANTS[0] / "jobs/sim1/generateIC.py",
                 "--config", job / "config.txt", "--seed", "42"])
        self.assertEqual(before, (job / "initial_conf.csv").read_bytes())
        data = rows(job / "initial_conf.csv")
        for axis in ("vx", "vy", "vz"):
            self.assertAlmostEqual(sum(float(p["m"]) * float(p[axis]) for p in data), 0, places=11)
        sources = [(ROOT / v / "jobs/sim1/generateIC.py").read_text() for v in VARIANTS]
        self.assertTrue(all(s == sources[0] for s in sources))

    def test_energy_and_snapshot_schema(self):
        for variant in VARIANTS:
            with self.subTest(variant=variant):
                job, cfg = self.job(variant)
                run_sim(variant, job, threads=4)
                energies = rows(job / "energy.csv")
                self.assertEqual([int(e["t"]) for e in energies], list(range(cfg["totsteps"]+1)))
                initial = rows(job / "OutputTimes/OutputTimes0.txt")
                final = rows(job / f"OutputTimes/OutputTimes{cfg['totsteps']}.txt")
                self.assertEqual(len(final), 64)
                self.assertTrue(all(len(p) == 18 and None not in p and None not in p.values() for p in final))
                for e in energies:
                    self.assertTrue(all(math.isfinite(float(x)) for x in e.values()))
                    self.assertAlmostEqual(float(e["K"])+float(e["U"]), float(e["Etot"]), places=10)
                box = [cfg[n] * cfg["rho"]**(-1/3) for n in ("nx", "ny", "nz")]
                for state, energy in ((initial, energies[0]), (final, energies[-1])):
                    independent = shifted_energy(state, box, cfg["rc"])
                    self.assertAlmostEqual(float(energy["U"]), independent, places=9)
                    self.assertAlmostEqual(sum(float(p["U"]) for p in state), independent, places=9)
                    for p in state:
                        for axis in "xyz":
                            self.assertAlmostEqual(float(p["a"+axis])*float(p["m"]), float(p["f"+axis]), places=10)
                drift = max(abs(float(e["Etot"])-float(energies[0]["Etot"])) for e in energies)
                self.assertLess(drift, 0.02 * abs(float(energies[0]["Etot"])))

    def test_serial_openmp_agreement(self):
        states = []
        for variant, threads in ((VARIANTS[1], 1), (VARIANTS[2], 1), (VARIANTS[2], 2), (VARIANTS[2], 4)):
            job, _ = self.job(f"{variant}_{threads}")
            run_sim(variant, job, threads)
            states.append(sorted(rows(job / "OutputTimes/OutputTimes20.txt"), key=lambda p: int(p["id"])))
        for state in states[1:]:
            for expected, actual in zip(states[0], state):
                for column in ("x", "y", "z", "vx", "vy", "vz", "fx", "fy", "fz", "U"):
                    self.assertAlmostEqual(float(actual[column]), float(expected[column]), places=10)

    def test_odd_small_and_reversed_mass_ratios_advance_full_time(self):
        for variant in VARIANTS:
            for mass1, mass2 in ((1.0, 9.0), (1.0, 2.25), (9.0, 1.0)):
                with self.subTest(variant=variant, mass1=mass1, mass2=mass2):
                    job, cfg = self.job(f"{variant}_{mass1}_{mass2}", nx=2, ny=1, nz=1,
                                        rho=0.002, rc=0.5, dt=0.1, totsteps=4,
                                        mass1=mass1, mass2=mass2, NPartRatio=0.5)
                    (job / "initial_conf.csv").write_text(
                        f"id,m,x,y,z,vx,vy,vz\n0,{min(mass1,mass2)},0,0,0,0.2,0.1,0\n"
                        f"1,{max(mass1,mass2)},3,0,0,-0.1,0,0.1\n")
                    run_sim(variant, job)
                    final = sorted(rows(job / "OutputTimes/OutputTimes4.txt"), key=lambda p: int(p["id"]))
                    self.assertAlmostEqual(float(final[0]["x"]), 0.08, places=12)
                    self.assertAlmostEqual(float(final[1]["x"]), 2.96, places=12)

    def test_time_reversibility(self):
        for variant in VARIANTS:
            for mass2 in (100.0, 9.0, 2.25):
                with self.subTest(variant=variant, mass2=mass2):
                    self.check_time_reversibility(variant, mass2)

    def check_time_reversibility(self, variant, mass2):
        job, cfg = self.job(f"{variant}_{mass2}", mass2=mass2)
        initial = sorted(rows(job / "initial_conf.csv"), key=lambda p: int(p["id"]))
        run_sim(variant, job)
        final = rows(job / "OutputTimes/OutputTimes20.txt")
        with (job / "initial_conf.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["id", "m", "x", "y", "z", "vx", "vy", "vz"])
            for p in final:
                writer.writerow([p[k] for k in ("id", "m", "x", "y", "z")] +
                                [-float(p[k]) for k in ("vx", "vy", "vz")])
        run_sim(variant, job)
        backward = sorted(rows(job / "OutputTimes/OutputTimes20.txt"), key=lambda p: int(p["id"]))
        for original, actual in zip(initial, backward):
            for key in ("x", "y", "z"):
                self.assertAlmostEqual(float(actual[key]), float(original[key]), places=10)
            for key in ("vx", "vy", "vz"):
                self.assertAlmostEqual(float(actual[key]), -float(original[key]), places=10)

    def test_tuckerman_disparate_mass_reference(self):
        # Exercise same-group and cross-group interactions over multiple macro
        # steps. An uncached full-force reference checks the optimized partial
        # buffers against the equations, rather than comparing the two engines.
        coordinates = ((0, 0, 0), (1.35, 0.1, 0),
                       (0.15, 1.4, 0.1), (1.4, 1.5, -0.1))
        velocities = ((0.05, 0.02, -0.01), (-0.03, 0.01, 0.02),
                      (0.01, -0.02, 0.03), (-0.02, -0.01, -0.02))
        for heavy_mass in (4.0, 100.0, 9.0, 2.25):
            initial = [[1.0 if i < 2 else heavy_mass, *coordinates[i], *velocities[i]]
                       for i in range(4)]
            cfg = dict(CONFIG, mass2=heavy_mass, NPartRatio=0.5,
                       nx=4, ny=1, nz=1, rho=0.001, dt=0.01, totsteps=7)
            expected_positions, expected_velocities = tuckerman_reference(initial, cfg)
            cases = (("MDMultipleTimeSteps", 1),
                     ("MDMultipleTimeStepsOpenMP", 1),
                     ("MDMultipleTimeStepsOpenMP", 4))
            for variant, threads in cases:
                with self.subTest(variant=variant, threads=threads, mass2=heavy_mass):
                    job = TEMP / self._testMethodName / f"{variant}_{threads}_{heavy_mass}"
                    job.mkdir(parents=True)
                    write_config(job / "config.txt", cfg)
                    with (job / "initial_conf.csv").open("w", newline="") as stream:
                        writer = csv.writer(stream)
                        writer.writerow(["id", "m", "x", "y", "z", "vx", "vy", "vz"])
                        for i, particle in enumerate(initial):
                            writer.writerow([i, *particle])
                    run_sim(variant, job, threads=threads)
                    actual = sorted(rows(job / "OutputTimes/OutputTimes7.txt"),
                                    key=lambda p: int(p["id"]))
                    self.assertEqual(len(actual), 4)
                    for i, particle in enumerate(actual):
                        for keys, expected in (("xyz", expected_positions),
                                               (("vx", "vy", "vz"), expected_velocities)):
                            for d, key in enumerate(keys):
                                self.assertAlmostEqual(float(particle[key]), expected[i][d],
                                                       delta=2e-12)

    def test_second_order_convergence(self):
        # Keep this pair inside rc throughout: cutoff crossings have a force
        # discontinuity and do not satisfy the smooth-force convergence premise.
        def smooth_job(name, dt, steps):
            job, _ = self.job(name, nx=2, ny=1, nz=1, rho=0.002,
                              dt=dt, totsteps=steps, NPartRatio=0.5)
            (job / "initial_conf.csv").write_text(
                "id,m,x,y,z,vx,vy,vz\n0,1,0,0,0,0.05,0.02,0\n"
                "1,100,1.3,0,0,-0.01,0,0\n")
            return job
        reference_job = smooth_job("reference", 0.0025, 80)
        run_sim(VARIANTS[0], reference_job)
        reference = sorted(rows(reference_job / "OutputTimes/OutputTimes80.txt"), key=lambda p: int(p["id"]))
        for variant in VARIANTS:
            errors = []
            for dt, steps in ((0.02, 10), (0.01, 20)):
                job = smooth_job(f"{variant}_{dt}", dt, steps)
                run_sim(variant, job)
                state = sorted(rows(job / f"OutputTimes/OutputTimes{steps}.txt"), key=lambda p: int(p["id"]))
                errors.append(math.sqrt(sum((float(p[k])-float(q[k]))**2
                                            for p, q in zip(state, reference)
                                            for k in ("x", "y", "z", "vx", "vy", "vz"))))
            with self.subTest(variant=variant):
                self.assertLess(errors[1], 0.4 * errors[0], f"Expected ~4x smaller error: {errors}")

    def test_zero_steps_and_single_species_normal(self):
        for variant in VARIANTS:
            job, _ = self.job(variant, totsteps=0)
            run_sim(variant, job)
            self.assertEqual(len(rows(job / "energy.csv")), 1)
        job, _ = self.job("one_species", mass1=1.0, mass2=1.0, totsteps=2)
        run_sim(VARIANTS[0], job)
        for variant in VARIANTS[1:]:
            run_sim(variant, job, succeeds=False)

    def test_invalid_configuration_fails_cleanly(self):
        invalid = ("dt -0.1", "dt nan", "dt 1oops", "mass1 0", "rho 0", "nx 0",
                   "nx 2.5", "NPartRatio 1.2", "totsteps -1", "rc 100", "mispelled 1", "dt")
        for variant in VARIANTS:
            for text in invalid:
                with self.subTest(variant=variant, invalid=text):
                    job, cfg = self.job(f"{variant}_{invalid.index(text)}")
                    key = text.split()[0]
                    contents = "\n".join(f"{k} {v}" for k, v in cfg.items() if k != key)
                    (job / "config.txt").write_text(contents + "\n" + text + "\n")
                    run_sim(variant, job, succeeds=False)
            job, cfg = self.job(f"{variant}_missing")
            del cfg["rho"]
            write_config(job / "config.txt", cfg)
            run_sim(variant, job, succeeds=False)
            job, _ = self.job(f"{variant}_duplicate")
            with (job / "config.txt").open("a") as stream:
                stream.write("dt 0.1\n")
            run_sim(variant, job, succeeds=False)

    def test_invalid_csv_fails_cleanly(self):
        for variant in VARIANTS:
            for kind in ("short", "nan", "count", "duplicate", "mass", "overlap", "header"):
                with self.subTest(variant=variant, kind=kind):
                    job, _ = self.job(f"{variant}_{kind}")
                    with (job / "initial_conf.csv").open() as stream:
                        data = list(csv.reader(stream))
                    if kind == "header": data[0][1] = "wrong"
                    if kind == "short": data[1] = ["0", "1"]
                    if kind == "nan": data[1][-1] = "nan"
                    if kind == "count": data.pop()
                    if kind == "duplicate": data[2][0] = data[1][0]
                    if kind == "mass": data[1][1] = "2.0"
                    if kind == "overlap": data[2][2:5] = data[1][2:5]
                    with (job / "initial_conf.csv").open("w", newline="") as stream:
                        csv.writer(stream).writerows(data)
                    run_sim(variant, job, succeeds=False)

    def test_mass_order_does_not_change_dynamics(self):
        for variant in VARIANTS:
            original, cfg = self.job(variant)
            run_sim(variant, original)
            expected = sorted(rows(original / "OutputTimes/OutputTimes20.txt"), key=lambda p: int(p["id"]))
            cfg["mass1"], cfg["mass2"] = cfg["mass2"], cfg["mass1"]
            write_config(original / "config.txt", cfg)
            run_sim(variant, original)
            actual = sorted(rows(original / "OutputTimes/OutputTimes20.txt"), key=lambda p: int(p["id"]))
            self.assertEqual(actual, expected)

    def test_nonfinite_final_step_fails(self):
        for variant in VARIANTS:
            with self.subTest(variant=variant):
                job, _ = self.job(variant, dt=1e308, totsteps=1)
                result = run_sim(variant, job, succeeds=False)
                self.assertIn("Non-finite", result.stderr)
                self.assertNotIn("All done", result.stdout)

    def test_missing_job_and_usage(self):
        for variant in VARIANTS:
            run_sim(variant, TEMP / "does_not_exist", succeeds=False)
            result = subprocess.run([str(BINARIES[variant])], capture_output=True, text=True)
            self.assertEqual(result.returncode, 1)
            self.assertIn("Usage:", result.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eigen-include", default=os.environ.get("EIGEN_INC", "/usr/include/eigen3"))
    parser.add_argument("--build-mode", choices=("debug", "releaseO2", "releaseO3"), default="debug")
    parser.add_argument("--sanitize", action="store_true", help="Compile temporary binaries with AddressSanitizer and UBSan")
    args = parser.parse_args()
    if not (Path(args.eigen_include) / "Eigen/Core").is_file():
        parser.error("Eigen/Core not found; set --eigen-include or EIGEN_INC")
    BINARIES, UNIT_BINARIES = {}, {}
    with tempfile.TemporaryDirectory(prefix="md-checks-") as directory:
        TEMP = Path(directory)
        for variant in VARIANTS:
            project = ROOT / variant
            flags = ["-std=c++17", "-Wall", "-Wextra", "-Wpedantic", "-I"+str(project/"src/include"), "-I"+args.eigen_include]
            if variant.endswith("OpenMP"): flags += ["-fopenmp"]
            cxx = shlex.split(os.environ.get("CXX", "g++"))
            if args.sanitize:
                flags += ["-O0", "-g", "-fsanitize=address,undefined", "-fno-omit-frame-pointer"]
                binary = TEMP / (variant + "_sanitized")
                # Compile once and share objects with the unit-check executable.
                # Keep sanitizer artifacts temporary, separate from normal builds.
                all_objects, objects = [], []
                for source in sorted((project/"src").rglob("*.cpp")):
                    obj = TEMP / variant / source.relative_to(project/"src").with_suffix(".o")
                    obj.parent.mkdir(parents=True, exist_ok=True)
                    command([*cxx, *flags, "-c", source, "-o", obj])
                    all_objects.append(obj)
                    if source.name != "MolecularDynamics.cpp": objects.append(obj)
                command([*cxx, *flags, *all_objects, "-o", binary, "-lstdc++fs"])
            else:
                command(["make", "-C", project, "-j2", "BUILD_MODE="+args.build_mode, "EIGEN_INC="+args.eigen_include])
                suffix = {"debug": "debug", "releaseO2": "O2", "releaseO3": "O3"}[args.build_mode]
                binary = project / "bin" / ("MolecularDynamics_"+suffix)
                objects = sorted(p for p in (project/"build"/args.build_mode).rglob("*.o") if p.name != "MolecularDynamics.o")
            BINARIES[variant] = binary
            unit = TEMP / (variant + "_unit")
            command([*cxx, *flags, ROOT/"validation/support/physics_checks.cpp", *objects, "-o", unit, "-lstdc++fs"])
            UNIT_BINARIES[variant] = unit
            print(f"Built {variant} ({'sanitized' if args.sanitize else args.build_mode})", flush=True)
        result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RegressionChecks))
        raise SystemExit(0 if result.wasSuccessful() else 1)
