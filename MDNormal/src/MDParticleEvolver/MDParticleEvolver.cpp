/**
 * @file MDParticleEvolver.cpp
 * @brief Implementation of the MDParticleEvolver class.
 *
 * Contains the core numerical engine: pairwise LJ force computation
 * using Newton's third law, Velocity Verlet time integration, and
 * state/energy output.
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

void MDParticleEvolver::evolve(double dt) {
    evolvePriv(dt);
}

double MDParticleEvolver::computeKineticEnergy(void) {
    return computeKineticEnergyPriv();
}

// ---------------------------------------------------------------------------
// Initialization
// ---------------------------------------------------------------------------

/**
 * @brief Identifies which mass is lighter and stores both values.
 */
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
 * Each particle is assigned the LJ parameters (sigma, epsilon, rc) and
 * box dimensions from the SimParams.  After loading, initial forces are
 * computed with a full pairwise sweep.
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

    U = 0;
    calcForcesPairwise();
}

// ---------------------------------------------------------------------------
// State output
// ---------------------------------------------------------------------------

/**
 * @brief Writes the full state (pos, vel, force, potential) of every
 *        particle to a CSV file.
 */
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

/**
 * @brief Appends a line of energy data (t, K, U, Etot) to a CSV log file.
 */
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
// Energy computation
// ---------------------------------------------------------------------------

/**
 * @brief Sums 0.5 * m * |v|^2 over all particles.
 */
double MDParticleEvolver::computeKineticEnergyPriv(void) {
    double kineticEnergy = 0;
    for (const auto& entry : particlesMap) {
        const auto& group = entry.second;
        for (const MDPart& particle : group) {
            kineticEnergy += 0.5 * particle.getMass()
                             * particle.getVelocity().squaredNorm();
        }
    }
    return kineticEnergy;
}

// ---------------------------------------------------------------------------
// Force computation
// ---------------------------------------------------------------------------

/**
 * @brief Computes all pairwise LJ forces using Newton's third law.
 *
 * 1. Flattens the mass-grouped map into a linear pointer array.
 * 2. Zeroes all forces.
 * 3. Iterates over all unique pairs (i, j) with j > i.
 * 4. Applies Newton's third law: F_ji = -F_ij.
 */
void MDParticleEvolver::calcForcesPairwise() {
    std::vector<MDPart*> allParticles;
    for (auto& entry : particlesMap) {
        for (MDPart& particle : entry.second) {
            particle.setForce(Eigen::Vector3d::Zero());
            allParticles.push_back(&particle);
        }
    }

    size_t numParticles = allParticles.size();

    for (size_t i = 0; i < numParticles; ++i) {
        for (size_t j = i + 1; j < numParticles; ++j) {
            Eigen::Vector3d force = allParticles[i]->computeForceijAndU(
                *allParticles[j], U);

            allParticles[i]->setForce(allParticles[i]->getForce() + force);
            allParticles[j]->setForce(allParticles[j]->getForce() - force);
        }
    }
}

// ---------------------------------------------------------------------------
// Time integration
// ---------------------------------------------------------------------------

/**
 * @brief Performs one full Velocity Verlet integration step.
 *
 * The Velocity Verlet algorithm is time-reversible and symplectic:
 *   1. v(t + dt/2) = v(t) + 0.5 * a(t) * dt
 *   2. r(t + dt)   = r(t) + v(t + dt/2) * dt
 *   3. Apply periodic boundary conditions
 *   4. Recompute forces → a(t + dt)
 *   5. v(t + dt)   = v(t + dt/2) + 0.5 * a(t + dt) * dt
 */
void MDParticleEvolver::VelocityVerlet(double dt) {
    // Half-kick + drift + PBC
    for (auto& entry : particlesMap) {
        auto& group = entry.second;
        for (MDPart& particle : group) {
            particle.updateVel(dt * 0.5);
            particle.updatePos(dt);
            particle.pbc();
        }
    }

    // Recompute forces at new positions
    U = 0;
    calcForcesPairwise();

    // Second half-kick
    for (auto& entry : particlesMap) {
        auto& group = entry.second;
        for (MDPart& particle : group) {
            particle.updateVel(dt * 0.5);
        }
    }
}

/**
 * @brief Evolves particles by one macro time-step using sub-stepping.
 *
 * The macro step dt is divided into Ratio sub-steps of dt/Ratio each.
 * This is the foundation of the multiple-time-step approach.
 */
void MDParticleEvolver::evolvePriv(double dt) {
    for (int i = 0; i < Ratio; i++) {
        VelocityVerlet(dt / Ratio);
    }
}
