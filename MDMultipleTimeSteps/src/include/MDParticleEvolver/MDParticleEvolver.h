/**
 * @file MDParticleEvolver.h
 * @brief Manages the particle collection and drives force computation
 *        and multiple-time-step (MTS) integration of particle states.
 *
 * MDParticleEvolver owns a map of particles grouped by mass and coordinates
 * the MTS Velocity Verlet integration scheme, pairwise force computation,
 * and I/O (state snapshots, energy logging).
 *
 * This is the serial (non-OpenMP) MTS variant.  For the parallelized
 * version, see the MDMultipleTimeStepsOpenMP variant.
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
 * @brief Orchestrates particle evolution with a serial multiple-time-step
 *        (MTS) Velocity Verlet integrator.
 *
 * In the MTS scheme, light particles are sub-stepped with a smaller dt
 * while heavy particles use the full dt.  This exploits the mass ratio
 * to reduce computational cost. The chosen scheme is Tuckerman, Berne and
 * Martyna (1992), Section IV.A, G_xyx with velocity-Verlet factors.
 * Both positions and velocities of the inactive group remain fixed.
 *
 * Same-mass contributions are saved in the auxiliary force buffer.
 * The primary force is reset, accumulated and combined as needed during
 * the MTS loop; its meaning depends on the current integration stage.
 */
class MDParticleEvolver {
public:
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

    /**
     * @brief Advances all particles by one macro time-step.
     *
     * Uses the MTS Velocity Verlet algorithm:
     *   1. Sub-step light particles for Ratio/2 sub-steps.
     *   2. Full step heavy particles.
     *   3. Sub-step light particles for Ratio/2 sub-steps again.
     * For odd Ratio, mirrored half small steps complete both light stages.
     */
    void evolve(std::vector<MDPart>& particleVecLight,
                std::vector<MDPart>& particleVecHeavy, double dt);

    void writeState(std::string path, std::string fileName, char delimiter);
    void writeEnergy(std::string path, std::string fileName,
                     long int t, char delimiter);
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
     * @brief Computes pairwise forces within a same-mass group.
     *
     * Uses Newton's third law (F_ji = -F_ij) and stores the result
     * in the auxiliary force buffer for later combination.
     */
    void calcForcesPairwiseSame(std::vector<MDPart>& particleVec);

    /**
     * @brief Computes pairwise forces between two different mass groups.
     */
    void calcForcesPairwiseDiff(std::vector<MDPart>& particleVec1,
                                std::vector<MDPart>& particleVec2);

    /**
     * @brief One Velocity Verlet step for a "free" particle group
     *        interacting with a "stuck" (stationary-for-this-step) group.
     */
    void VelocityVerlet(std::vector<MDPart>& particleVecFree,
                        std::vector<MDPart>& particleVecStuck, double dt);

    /**
     * @brief Combines primary and auxiliary forces for a particle group.
     *
     * Sets total force = same-mass force + cross-mass force.
     */
    void combineForces(std::vector<MDPart>& particleVec);

    /// Private implementation of the MTS evolution step.
    void evolvePriv(std::vector<MDPart>& particleVecLight,
                    std::vector<MDPart>& particleVecHeavy, double dt);
};

#endif /* MDPARTICLEEVOLVER_H_ */
