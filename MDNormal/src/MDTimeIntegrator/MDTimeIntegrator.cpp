/**
 * @file MDTimeIntegrator.cpp
 * @brief Implementation of the MDTimeIntegrator class.
 *
 * Contains the main simulation loop that drives the Molecular Dynamics
 * simulation from t = 0 to t = totsteps * dt.
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
 * @brief Main simulation loop.
 *
 * For each time step t in [0, totsteps]:
 *   - Every 100 steps and at the final step: write a full state snapshot to a file.
 *   - Every step: append energy data to the energy CSV.
 *   - Unless this is the final state, evolve particles by dt using the Velocity Verlet integrator.
 *
 * The entire loop is wrapped in a try-catch to handle runtime errors
 * gracefully (file I/O failures and explicit non-finite-state checks).
 */
bool MDTimeIntegrator::runPriv(void) {
    try {
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

            // Advance the simulation
            this->evolver->evolve(dt);
        }

        return true;

    } catch (const std::exception& e) {
        std::cerr << "Simulation crashed: " << e.what() << std::endl;
        return false;
    }
}
