/**
 * @file MDTimeIntegrator.h
 * @brief Drives the main simulation loop, advancing particles through
 *        time and coordinating I/O output at regular intervals.
 *
 * @author Federico Celauro
 * @date   Apr 6, 2026
 */

#ifndef INCLUDE_MDTIMEINTEGRATOR_MDTIMEINTEGRATOR_H_
#define INCLUDE_MDTIMEINTEGRATOR_MDTIMEINTEGRATOR_H_

#include "../MDParticleEvolver/MDParticleEvolver.h"

/**
 * @class MDTimeIntegrator
 * @brief Controls the simulation time loop and delegates particle
 *        evolution to MDParticleEvolver.
 *
 * At each time step the integrator:
 *   1. Optionally writes a full state snapshot (every 100 steps).
 *   2. Logs energy (kinetic, potential, total).
 *   3. Advances the particle evolver by dt.
 */
class MDTimeIntegrator {
public:
    MDTimeIntegrator(SimParams params, MDParticleEvolver* evolver);
    virtual ~MDTimeIntegrator();

    /// @name Time-step accessors
    /// @{
    long int getTotsteps() const;
    void setTotsteps(long int totsteps);
    double getDt() const;
    void setDt(double dt);
    /// @}

    /// @name Evolver accessors
    /// @{
    MDParticleEvolver*& getEvolver();
    void setEvolver(MDParticleEvolver*& evolver);
    /// @}

    /// @name Output path accessors
    /// @{
    const std::string& getOutputTimePath() const;
    void setOutputTimePath(const std::string& path);
    const std::string& getEnergyPath() const;
    void setEnergyPath(const std::string& path);
    /// @}

    /**
     * @brief Executes the full simulation loop.
     * @return true if the simulation completed successfully, false otherwise.
     */
    bool run(void);

private:
    long int totsteps;            ///< Total number of integration steps
    double dt;                    ///< Time-step size
    std::string OutputTimePath;   ///< Directory for state snapshot files
    std::string EnergyPath;       ///< Directory for the energy log file
    MDParticleEvolver* evolver;   ///< Non-owning pointer to the particle evolver

    /// Private implementation of the simulation loop.
    bool runPriv(void);
};

#endif /* INCLUDE_MDTIMEINTEGRATOR_MDTIMEINTEGRATOR_H_ */
