/**
 * @file MDParticleEvolver.cpp
 * @brief Implementation of the MDParticleEvolver class (OpenMP variant).
 *
 * This is the computational core of the OpenMP-parallelized
 * multiple-time-step Molecular Dynamics simulation.  Key parallelized
 * sections include:
 *   - Pairwise force computation (parallel for + atomic accumulation)
 *   - Kinetic energy summation (parallel for + reduction)
 *   - Velocity Verlet kicks and drifts (parallel for)
 *
 * Reproduces the disparate-mass G_xyx splitting of Tuckerman, Berne and
 * Martyna (1992), Section IV.A, Eqs. (4.2)-(4.6); DOI: 10.1063/1.463137.
 * The inactive group's positions AND velocities are intentionally fixed.
 * See docs/ALGORITHM.md and docs/IMPLEMENTATION.md in the repository root.
 * This coordinate-group splitting is reversible and second order for smooth
 * forces, but is not generally canonical symplectic for coupled mass groups.
 * See docs/ALGORITHM.md for the energy-error comparison and method properties.
 *
 * @author Federico Celauro
 * @date   Apr 6, 2026
 */

#include "../include/MDParticleEvolver/MDParticleEvolver.h"
#include <set>
#include <cmath>
#include <iomanip>
#include <limits>
#include <algorithm>

// ---------------------------------------------------------------------------
// Constructor & Destructor
// ---------------------------------------------------------------------------

MDParticleEvolver::MDParticleEvolver(SimParams params, std::string path,
                                     std::string fileName, char delimiter) {
    initializePriv(params, path, fileName, delimiter);
    initializeMass(params.getMass1(), params.getMass2());
    setRatio(params.getRatio());
}

MDParticleEvolver::~MDParticleEvolver() {
}

// ---------------------------------------------------------------------------
// Mass accessors
// ---------------------------------------------------------------------------

double MDParticleEvolver::getHeavyMass() const { return heavyMass; }
void MDParticleEvolver::setHeavyMass(double heavyMass) { this->heavyMass = heavyMass; }

double MDParticleEvolver::getLightMass() const { return lightMass; }
void MDParticleEvolver::setLightMass(double lightMass) { this->lightMass = lightMass; }

int MDParticleEvolver::getRatio() const { return Ratio; }
void MDParticleEvolver::setRatio(int ratio) { Ratio = ratio; }

// ---------------------------------------------------------------------------
// Potential energy accessors
// ---------------------------------------------------------------------------

double MDParticleEvolver::getU() const { return U; }
void MDParticleEvolver::setU(double u) { U = u; }

double MDParticleEvolver::getUAux() const { return UAux; }
void MDParticleEvolver::setUAux(double uAux) { UAux = uAux; }

// ---------------------------------------------------------------------------
// Particle map accessors
// ---------------------------------------------------------------------------

std::map<double, std::vector<MDPart>>& MDParticleEvolver::getParticlesMap() {
    return particlesMap;
}

const std::map<double, std::vector<MDPart>>& MDParticleEvolver::getParticlesMap() const {
    return particlesMap;
}

void MDParticleEvolver::setParticlesMap(
        const std::map<double, std::vector<MDPart>>& particlesMap) {
    this->particlesMap = particlesMap;
}

// ---------------------------------------------------------------------------
// Public simulation interface
// ---------------------------------------------------------------------------

void MDParticleEvolver::writeState(std::string path, std::string fileName, char delimiter) {
    writeStatePriv(path, fileName, delimiter);
}

void MDParticleEvolver::writeEnergy(std::string path, std::string fileName,
                                     long int t, char delimiter) {
    writeEnergyPriv(path, fileName, t, delimiter);
}

void MDParticleEvolver::evolve(std::vector<MDPart>& particleVecLight,
                                std::vector<MDPart>& particleVecHeavy, double dt) {
    evolvePriv(particleVecLight, particleVecHeavy, dt);
}

double MDParticleEvolver::computeKineticEnergy(void) {
    return computeKineticEnergyPriv();
}

// ---------------------------------------------------------------------------
// Initialization
// ---------------------------------------------------------------------------

void MDParticleEvolver::initializeMass(double mass1, double mass2) {
    if (mass1 <= mass2) {
        setLightMass(mass1);
        setHeavyMass(mass2);
    } else {
        setLightMass(mass2);
        setHeavyMass(mass1);
    }
}

/**
 * @brief Reads initial conditions from a CSV file and groups particles by mass.
 *
 * After loading, initial forces are calculated with a full pairwise sweep
 * (same-mass + cross-mass) and combined for the light particles.
 */
void MDParticleEvolver::initializePriv(SimParams params, std::string path,
                                        std::string fileName, char delimiter = ',') {
    asciiFileReader fileReader(path, fileName);
    std::vector<std::string> csvLines = fileReader.getAllLines();
    if (csvLines.empty() ||
        (csvLines.front() != "id,m,x,y,z,vx,vy,vz" &&
         csvLines.front() != "id,m,x,y,z,vx,vy,vz\r")) {
        throw std::invalid_argument("Expected CSV header: id,m,x,y,z,vx,vy,vz");
    }
    csvLines.erase(csvLines.begin());

    std::set<long long int> ids;
    for (const std::string& line : csvLines) {
        if (line.find_first_not_of(" \t\r") == std::string::npos) continue;

        MDPart newPart(line, delimiter);
        double mass = newPart.getMass();
        if (mass != params.getMass1() && mass != params.getMass2()) {
            throw std::invalid_argument("CSV mass does not match either configured mass.");
        }
        if (!ids.insert(newPart.getId()).second) {
            throw std::invalid_argument("Duplicate particle ID in CSV.");
        }

        newPart.setRc(params.getRc());
        newPart.setSigma(params.getSigma());
        newPart.setEpsilon(params.getEpsilon());
        newPart.setLx(params.getLx());
        newPart.setLy(params.getLy());
        newPart.setLz(params.getLz());
        newPart.setUcut();

        newPart.pbc();
        particlesMap[mass].push_back(newPart);
    }

    if (ids.size() != static_cast<size_t>(params.getNparticles())) {
        throw std::invalid_argument("CSV particle count must equal nx * ny * nz. Regenerate initial conditions after changing config.txt.");
    }
    // Reject periodic overlaps before entering an OpenMP region, where throwing is unsafe.
    std::vector<const MDPart*> loaded;
    for (const auto& entry : particlesMap) {
        for (const auto& particle : entry.second) loaded.push_back(&particle);
    }
    const Eigen::Vector3d box(params.getLx(), params.getLy(), params.getLz());
    for (size_t i = 0; i < loaded.size(); ++i) {
        for (size_t j = i + 1; j < loaded.size(); ++j) {
            Eigen::Vector3d rij = loaded[i]->getPosition() - loaded[j]->getPosition();
            rij.array() -= box.array() * (rij.array() / box.array()).round();
            if (rij.squaredNorm() == 0) throw std::invalid_argument("Overlapping particles in initial conditions.");
        }
    }

    if (params.getMass1() == params.getMass2() || particlesMap.size() != 2) {
        throw std::invalid_argument("MTS requires two distinct, nonempty mass groups; use MDNormal for one species.");
    }
    std::vector<MDPart>& LightPartVec = getParticlesMap().at(std::min(params.getMass1(), params.getMass2()));
    std::vector<MDPart>& HeavyPartVec = getParticlesMap().at(std::max(params.getMass1(), params.getMass2()));

    // Compute initial forces for each mass group and cross-group interactions
    calcForcesPairwiseSame(LightPartVec);
    calcForcesPairwiseSame(HeavyPartVec);
    calcForcesPairwiseDiff(HeavyPartVec, LightPartVec);
    combineForces(LightPartVec);
}

// ---------------------------------------------------------------------------
// State output
// ---------------------------------------------------------------------------

void MDParticleEvolver::writeStatePriv(std::string path, std::string fileName,
                                        char delimiter) {
    asciiFileWriter fileWriter(path, fileName);
    std::string header = std::string("id") + delimiter + "m" + delimiter
        + "x" + delimiter + "y" + delimiter + "z" + delimiter
        + "vx" + delimiter + "vy" + delimiter + "vz" + delimiter
        + "ax" + delimiter + "ay" + delimiter + "az" + delimiter
        + "fx" + delimiter + "fy" + delimiter + "fz" + delimiter
        + "fxAux" + delimiter + "fyAux" + delimiter + "fzAux" + delimiter + "U";

    fileWriter.writeOnFile(header);
    fileWriter.setAppendMode();

    // Recompute diagnostics on copies: MTS buffers contain partial forces on
    // the stationary group and must not be overwritten by snapshot output.
    std::vector<MDPart> snapshot;
    for (const auto& entry : particlesMap) {
        for (const auto& particle : entry.second) {
            snapshot.push_back(particle);
            snapshot.back().setForce(Eigen::Vector3d::Zero());
            snapshot.back().setU(0.0);
        }
    }
    for (size_t i = 0; i < snapshot.size(); ++i) {
        for (size_t j = i + 1; j < snapshot.size(); ++j) {
            double energy = 0.0;
            const Eigen::Vector3d force = snapshot[i].computeForceijAndU(snapshot[j], energy);
            snapshot[i].setForce(snapshot[i].getForce() + force);
            snapshot[j].setForce(snapshot[j].getForce() - force);
            snapshot[i].setU(snapshot[i].getU() + energy * 0.5);
            snapshot[j].setU(snapshot[j].getU() + energy * 0.5);
        }
    }
    for (const auto& particle : snapshot) {
        fileWriter.writeOnFile(particle.returnStateString(delimiter));
    }
}

void MDParticleEvolver::writeEnergyPriv(std::string path, std::string fileName,
                                         long int t, char delimiter) {
    double K;
    asciiFileWriter fileWriter(path, fileName, true);

    K = computeKineticEnergy();
    if (!std::isfinite(K) || !std::isfinite(getU())) {
        throw std::runtime_error("Non-finite energy; reduce dt and check the initial conditions.");
    }
    for (const auto& entry : particlesMap) {
        for (const auto& particle : entry.second) {
            if (!particle.getPosition().allFinite() || !particle.getVelocity().allFinite()) {
                throw std::runtime_error("Non-finite particle state; reduce dt.");
            }
        }
    }
    fileWriter.getPointer() << std::setprecision(std::numeric_limits<double>::max_digits10)
                            << t << delimiter
                            << K << delimiter
                            << getU() << delimiter
                            << K + getU() << std::endl;
}

// ---------------------------------------------------------------------------
// Kinetic energy computation (OpenMP reduction)
// ---------------------------------------------------------------------------

/**
 * @brief Computes total kinetic energy using OpenMP reduction.
 *
 * The outer parallel region distributes over mass groups, and the
 * inner `#pragma omp for nowait` distributes over particles within
 * each group.  The `nowait` clause allows threads to proceed without
 * a barrier between groups.
 */
double MDParticleEvolver::computeKineticEnergyPriv(void) {
    double kineticEnergy = 0;

    #pragma omp parallel reduction(+:kineticEnergy)
    for (const auto& entry : particlesMap) {
        const auto& group = entry.second;

        #pragma omp for nowait
        for (size_t i = 0; i < group.size(); ++i) {
            const MDPart& particle = group[i];
            kineticEnergy += 0.5 * particle.getMass()
                             * particle.getVelocity().squaredNorm();
        }
    }
    return kineticEnergy;
}

// ---------------------------------------------------------------------------
// Pairwise force computation (OpenMP parallelized)
// ---------------------------------------------------------------------------

/**
 * @brief Computes all pairwise forces within a single mass group.
 *
 * Parallelization strategy:
 *   1. Zero all forces in parallel.
 *   2. Iterate over unique pairs (i, j) with j > i using `#pragma omp parallel for`.
 *   3. Use addForceAtomic() (omp atomic) for thread-safe force updates.
 *   4. Use `reduction(+:totalU)` for the potential energy sum.
 *   5. Save the same-group force as the auxiliary force for later combination.
 */
void MDParticleEvolver::calcForcesPairwiseSame(std::vector<MDPart>& particleVec) {
    // Zero all forces
    #pragma omp parallel for
    for (size_t i = 0; i < particleVec.size(); ++i) {
        particleVec[i].setForce(Eigen::Vector3d::Zero());
    }

    double totalU = 0;

    // Compute pairwise interactions
    #pragma omp parallel for reduction(+:totalU)
    for (size_t i = 0; i < particleVec.size(); ++i) {
        MDPart& particle1 = particleVec[i];

        for (size_t j = i + 1; j < particleVec.size(); ++j) {
            MDPart& particle2 = particleVec[j];

            double energy = 0;
            Eigen::Vector3d force = particle1.computeForceijAndU(particle2, energy);

            totalU += energy;
            particle1.addForceAtomic(force);
            particle2.addForceAtomic(-force);
        }
    }
    this->U += totalU;

    // Store same-group force as auxiliary for later MTS combination
    #pragma omp parallel for
    for (size_t i = 0; i < particleVec.size(); ++i) {
        particleVec[i].setForceAux(particleVec[i].getForce());
        particleVec[i].setAccelerationFromForce();
    }
}

/**
 * @brief Computes pairwise forces between two different mass groups.
 *
 * particleVec1 is the "outer" group (parallelized over), and
 * particleVec2 is the "inner" group.  Forces are accumulated
 * atomically on both groups using Newton's third law.
 */
void MDParticleEvolver::calcForcesPairwiseDiff(std::vector<MDPart>& particleVec1,
                                                std::vector<MDPart>& particleVec2) {
    // Zero forces on the second group (first group already zeroed in calcForcesPairwiseSame)
    #pragma omp parallel for
    for (size_t i = 0; i < particleVec2.size(); ++i) {
        particleVec2[i].setForce(Eigen::Vector3d::Zero());
    }

    double totalU = 0;

    #pragma omp parallel for reduction(+:totalU)
    for (size_t i = 0; i < particleVec1.size(); ++i) {
        MDPart& particle1 = particleVec1[i];

        for (size_t j = 0; j < particleVec2.size(); ++j) {
            MDPart& particle2 = particleVec2[j];

            double energy = 0;
            Eigen::Vector3d force = particle1.computeForceijAndU(particle2, energy);

            totalU += energy;
            particle1.addForceAtomic(force);
            particle2.addForceAtomic(-force);
        }
    }
    this->U += totalU;

    #pragma omp parallel for
    for (size_t i = 0; i < particleVec1.size(); ++i) {
        particleVec1[i].setAccelerationFromForce();
    }
    #pragma omp parallel for
    for (size_t i = 0; i < particleVec2.size(); ++i) {
        particleVec2[i].setAccelerationFromForce();
    }
}

/**
 * @brief Combines primary and auxiliary forces.
 *
 * The total force on a particle is the sum of the same-group force
 * (fAux, from calcForcesPairwiseSame) and the cross-group force
 * (f, from calcForcesPairwiseDiff).
 */
void MDParticleEvolver::combineForces(std::vector<MDPart>& particleVec) {
    #pragma omp parallel for
    for (size_t i = 0; i < particleVec.size(); ++i) {
        MDPart& particle = particleVec[i];
        particle.setForce(particle.getForce() + particle.getForceAux());
    }
}

// ---------------------------------------------------------------------------
// Velocity Verlet integration (OpenMP parallelized)
// ---------------------------------------------------------------------------

/**
 * @brief Performs one Velocity Verlet step for the "free" particle group.
 *
 * The "free" particles are updated while the "stuck" group remains
 * fixed in both position and velocity for this sub-step. Opposite cross
 * forces are cached on that group, but no inactive momentum kick is applied;
 * this is the selected paper scheme, not a Hamiltonian potential split.
 *
 * Steps (all parallelized with OpenMP):
 *   1. Half-kick:  v += 0.5 * a * dt
 *   2. Drift:      r += v * dt, then apply PBC
 *   3. Recompute forces (same-mass + cross-mass)
 *   4. Half-kick:  v += 0.5 * a * dt
 */
void MDParticleEvolver::VelocityVerlet(std::vector<MDPart>& particleVecFree,
                                        std::vector<MDPart>& particleVecStuck,
                                        double dt) {
    // Step 1-2: Half-kick + drift + PBC
    #pragma omp parallel for
    for (size_t i = 0; i < particleVecFree.size(); ++i) {
        MDPart& particle = particleVecFree[i];
        particle.updateVel(dt * 0.5);
        particle.updatePos(dt);
        particle.pbc();
    }

    // Step 3: Recompute forces at new positions
    U = 0;
    calcForcesPairwiseSame(particleVecFree);
    if (particleVecFree[0].getMass() == heavyMass) UAux = U;

    calcForcesPairwiseDiff(particleVecFree, particleVecStuck);

    // Step 4: Second half-kick
    #pragma omp parallel for
    for (size_t i = 0; i < particleVecFree.size(); ++i) {
        particleVecFree[i].updateVel(dt * 0.5);
    }
}

// ---------------------------------------------------------------------------
// Multiple-time-step evolution
// ---------------------------------------------------------------------------

/**
 * @brief Multiple-time-step (MTS) Velocity Verlet integration.
 *
 * The algorithm is symmetric about the heavy-particle step:
 *   1. Sub-step light particles for Ratio/2 sub-steps (dt_light = dt/Ratio)
 *      For an odd Ratio, add a half small step on each side of the heavy step.
 *   2. Combine forces on heavy particles
 *   3. Full Velocity Verlet step for heavy particles (dt_heavy = dt)
 *   4. Combine forces on light particles
 *   5. Sub-step light particles for another Ratio/2 sub-steps
 *   6. Add auxiliary potential energy to total
 *
 * This ensures that light particles (which move faster) see updated
 * forces more frequently, while heavy particles use a larger time step.
 */
void MDParticleEvolver::evolvePriv(std::vector<MDPart>& particleVecLight,
                                    std::vector<MDPart>& particleVecHeavy,
                                    double dt) {
    // First half of light-particle sub-steps
    for (int i = 0; i < Ratio / 2; i++) {
        VelocityVerlet(particleVecLight, particleVecHeavy, dt / Ratio);
    }

    if (Ratio % 2 != 0) {
        VelocityVerlet(particleVecLight, particleVecHeavy, dt / (2.0 * Ratio));
    }

    // Heavy-particle full step
    UAux = 0;
    combineForces(particleVecHeavy);
    VelocityVerlet(particleVecHeavy, particleVecLight, dt);
    combineForces(particleVecLight);

    // Reverse the first-half step order to keep odd ratios time-symmetric.
    if (Ratio % 2 != 0) {
        VelocityVerlet(particleVecLight, particleVecHeavy, dt / (2.0 * Ratio));
    }
    for (int i = 0; i < Ratio / 2; i++) {
        VelocityVerlet(particleVecLight, particleVecHeavy, dt / Ratio);
    }

    U += UAux;
}
