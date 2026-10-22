/**
 * @file MDTimeIntegrator.cpp
 * @brief Implementation of the MDTimeIntegrator class (MTS variant).
 *
 * Contains the main simulation loop that drives the multiple-time-step
 * Molecular Dynamics simulation.  Particle vectors for light and heavy
 * species are passed to the evolver at each step.
 *
 * @author Federico Celauro
 * @date   Apr 6, 2026
 */

#include "../include/MDTimeIntegrator/MDTimeIntegrator.h"

// ---------------------------------------------------------------------------
// Constructor & Destructor
// ---------------------------------------------------------------------------

MDTimeIntegrator::MDTimeIntegrator(SimParams params, MDParticleEvolver* evolver) {
    setEvolver(evolver);
    setDt(params.getDt());
    setTotsteps(params.getTotsteps());
    setOutputTimePath(params.getOutputTimePath());
    setEnergyPath(params.getEnergyPath());
}

MDTimeIntegrator::~MDTimeIntegrator() {
}

// ---------------------------------------------------------------------------
// Time-step accessors
// ---------------------------------------------------------------------------

long int MDTimeIntegrator::getTotsteps() const { return totsteps; }
void MDTimeIntegrator::setTotsteps(long int totsteps) { this->totsteps = totsteps; }

double MDTimeIntegrator::getDt() const { return dt; }
void MDTimeIntegrator::setDt(double dt) { this->dt = dt; }

// ---------------------------------------------------------------------------
// Evolver accessors
// ---------------------------------------------------------------------------

MDParticleEvolver*& MDTimeIntegrator::getEvolver() { return evolver; }
void MDTimeIntegrator::setEvolver(MDParticleEvolver*& evolver) { this->evolver = evolver; }

// ---------------------------------------------------------------------------
// Output path accessors
// ---------------------------------------------------------------------------

const std::string& MDTimeIntegrator::getOutputTimePath() const { return OutputTimePath; }
void MDTimeIntegrator::setOutputTimePath(const std::string& outputTimePath) {
    OutputTimePath = outputTimePath;
}

const std::string& MDTimeIntegrator::getEnergyPath() const { return EnergyPath; }
void MDTimeIntegrator::setEnergyPath(const std::string& energyPath) {
    EnergyPath = energyPath;
}

// ---------------------------------------------------------------------------
// Simulation loop
// ---------------------------------------------------------------------------

bool MDTimeIntegrator::run(void) {
    return runPriv();
}

/**
 * @brief Main simulation loop (MTS variant).
 *
 * Retrieves the light and heavy particle vectors from the evolver and
 * passes them explicitly to evolve() at each step.
 *
 * For each time step t in [0, totsteps]:
 *   - Every 100 steps and at the final step: write a full state snapshot.
 *   - Every step: append energy data to the energy CSV.
 *   - Unless this is the final state, evolve particles with the MTS Velocity Verlet integrator.
 */
bool MDTimeIntegrator::runPriv(void) {
    try {
        double lightMass = this->evolver->getLightMass();
        double heavyMass = this->evolver->getHeavyMass();

        // Validate that both mass groups exist in the particle map
        std::map<double, std::vector<MDPart>>& particlesMap =
            this->evolver->getParticlesMap();

        if (particlesMap.find(lightMass) == particlesMap.end() ||
            particlesMap.find(heavyMass) == particlesMap.end()) {
            std::cerr << "Error: Mass keys not found in particles map."
                      << std::endl;
            return false;
        }

        std::vector<MDPart>& LightPartVec =
            this->evolver->getParticlesMap().at(lightMass);
        std::vector<MDPart>& HeavyPartVec =
            this->evolver->getParticlesMap().at(heavyMass);

        // Main time-stepping loop
        for (long int t = 0; t <= totsteps; t++) {
            // Write state snapshot every 100 steps
            if (t % 100 == 0 || t == totsteps) {
                std::string fileName = "OutputTimes" + std::to_string(t) + ".txt";
                this->evolver->writeState(OutputTimePath, fileName, ',');
            }

            // Log energy at every step
            this->evolver->writeEnergy(EnergyPath, "energy.csv", t, ',');

            if (t == totsteps) break;

            // Advance the simulation via MTS Velocity Verlet
            this->evolver->evolve(LightPartVec, HeavyPartVec, dt);
        }

        return true;

    } catch (const std::exception& e) {
        std::cerr << "Simulation crashed: " << e.what() << std::endl;
        return false;
    }
}
