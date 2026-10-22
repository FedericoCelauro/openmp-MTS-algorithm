/**
 * @file MolecularDynamics.h
 * @brief Top-level class that initializes and runs a Molecular Dynamics
 *        simulation.
 *
 * MolecularDynamics ties together the simulation parameters, particle
 * evolver, and time integrator into a single entry point.
 *
 * @author Federico Celauro
 * @date   Apr 8, 2026
 */

#ifndef INCLUDE_MOLECULARDYNAMICS_H_
#define INCLUDE_MOLECULARDYNAMICS_H_

#include <iostream>
#include <string>
#include <memory>
#include <chrono>

#include "MDParticleEvolver/MDParticleEvolver.h"
#include "MDTimeIntegrator/MDTimeIntegrator.h"

/**
 * @class MolecularDynamics
 * @brief Orchestrates the full MD simulation pipeline.
 *
 * Workflow:
 *   1. Read simulation parameters from a config file.
 *   2. Initialize particles from an initial-conditions CSV.
 *   3. Create the time integrator and run the simulation loop.
 */
class MolecularDynamics {
public:
    MolecularDynamics();
    MolecularDynamics(char** argv);
    virtual ~MolecularDynamics();

    /// @name Component accessors
    /// @{
    MDParticleEvolver* getParticleEvolver() const;
    void setParticleEvolver(std::unique_ptr<MDParticleEvolver> evolver);

    MDTimeIntegrator* getTimeIntegrator() const;
    void setTimeIntegrator(std::unique_ptr<MDTimeIntegrator> integrator);
    /// @}

    /// @name Path accessors
    /// @{
    const std::string& getPathInitialConditions() const;
    void setPathInitialConditions(const std::string& path);

    const std::string& getConfigsFileName() const;
    void setConfigsFileName(const std::string& configsFileName);

    const std::string& getInitialConditionsFileName() const;
    void setInitialConditionsFileName(const std::string& fileName);
    /// @}

    /// @name Simulation name accessors
    /// @{
    const std::string& getSimName() const;
    void setSimName(const std::string& simName);
    /// @}

    /**
     * @brief Starts the simulation.
     * @return true if the simulation completed successfully.
     */
    bool start();

private:
    std::string SimName;
    std::string pathInitialConditions;
    std::string initialConditionsFileName = "initial_conf.csv";
    std::string configsFileName = "config.txt";

    /// Particle evolver (owned)
    std::unique_ptr<MDParticleEvolver> particleEvolver;

    /// Time integrator (owned)
    std::unique_ptr<MDTimeIntegrator> timeIntegrator;

    /**
     * @brief Initializes all simulation components.
     *
     * Reads parameters, creates the particle evolver and time
     * integrator, and writes the energy CSV header.
     *
     * @return true if initialization was successful.
     */
    bool init();
};

#endif /* INCLUDE_MOLECULARDYNAMICS_H_ */
