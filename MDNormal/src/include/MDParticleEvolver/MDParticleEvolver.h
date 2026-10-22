/**
 * @file MDParticleEvolver.h
 * @brief Manages the particle collection and drives force computation
 *        and time integration of particle states.
 *
 * MDParticleEvolver owns a map of particles grouped by mass and coordinates
 * the Velocity Verlet integration scheme, pairwise force computation,
 * and I/O (state snapshots, energy logging).
 *
 * @author Federico Celauro
 * @date   Apr 6, 2026
 */

#ifndef MDPARTICLEEVOLVER_H_
#define MDPARTICLEEVOLVER_H_

#include <string>
#include <map>
#include "../IO/ASCII/asciiFileWriter.h"
#include "../MDPart/MDPart.h"
#include "../SimParams/SimParams.h"

/**
 * @class MDParticleEvolver
 * @brief Orchestrates particle evolution: force computation, integration,
 *        and output writing.
 *
 * Particles are stored in a std::map keyed by mass, allowing the
 * multiple-time-step variants to iterate over mass groups independently.
 */
class MDParticleEvolver {
public:
    /**
     * @brief Constructs the evolver by reading initial conditions from file.
     * @param params     Simulation parameters.
     * @param path       Directory of the initial conditions file.
     * @param fileName   Name of the CSV initial conditions file.
     * @param delimiter  CSV delimiter character.
     */
    MDParticleEvolver(SimParams params, std::string path,
                      std::string fileName, char delimiter);

    virtual ~MDParticleEvolver();

    /// @name Mass accessors
    /// @{
    double getHeavyMass() const;
    void setHeavyMass(double heavyMass);
    double getLightMass() const;
    void setLightMass(double lightMass);
    int getRatio() const;
    void setRatio(int ratio);
    /// @}

    /// @name Potential energy accessors
    /// @{
    double getU() const;
    void setU(double u);
    double getUAux() const;
    void setUAux(double uAux);
    /// @}

    /// @name Particle map accessors
    /// @{
    std::map<double, std::vector<MDPart>>& getParticlesMap();
    const std::map<double, std::vector<MDPart>>& getParticlesMap() const;
    void setParticlesMap(const std::map<double, std::vector<MDPart>>& particlesMap);
    /// @}

    /// @name Simulation operations
    /// @{

    /// Advances all particles by one macro time-step dt.
    void evolve(double dt);

    /// Writes the full particle state to a CSV file.
    void writeState(std::string path, std::string fileName, char delimiter);

    /// Appends kinetic, potential, and total energy to a CSV log.
    void writeEnergy(std::string path, std::string fileName,
                     long int t, char delimiter);

    /// Computes the total kinetic energy of all particles.
    double computeKineticEnergy(void);
    /// @}

private:
    double lightMass, heavyMass;  ///< Light and heavy particle masses
    double U = 0.0, UAux = 0.0;               ///< Total and auxiliary potential energy
    int Ratio;                    ///< Time-step ratio sqrt(m_heavy / m_light)

    /// Particles grouped by mass
    std::map<double, std::vector<MDPart>> particlesMap;

    /// @name Initialization
    /// @{
    void initializeMass(double mass1, double mass2);
    void initializePriv(SimParams params, std::string path,
                        std::string fileName, char delimiter);
    /// @}

    /// @name Private implementations
    /// @{
    void writeStatePriv(std::string path, std::string fileName, char delimiter);
    void writeEnergyPriv(std::string path, std::string fileName,
                         long int t, char delimiter);
    double computeKineticEnergyPriv(void);
    /// @}

    /**
     * @brief Computes all pairwise LJ forces using Newton's third law.
     *
     * Flattens the particle map into a linear array, then iterates over
     * all unique pairs (i, j) with j > i.
     */
    void calcForcesPairwise();

    /**
     * @brief One full Velocity Verlet step.
     *
     * 1. Half-kick:  v += 0.5 * a * dt
     * 2. Drift:      r += v * dt, then apply PBC
     * 3. Recompute forces
     * 4. Half-kick:  v += 0.5 * a * dt
     */
    void VelocityVerlet(double dt);

    /// Private implementation of evolve().
    void evolvePriv(double dt);
};

#endif /* MDPARTICLEEVOLVER_H_ */
