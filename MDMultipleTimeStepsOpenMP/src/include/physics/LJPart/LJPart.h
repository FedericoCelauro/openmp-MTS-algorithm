/**
 * @file LJPart.h
 * @brief Extends MassivePart with Lennard-Jones interaction forces
 *        and thread-safe atomic force accumulation (OpenMP variant).
 *
 * An LJPart is a massive particle that interacts with other particles
 * via the Lennard-Jones 12-6 potential:
 *   U(r) = 4 * epsilon * [ (sigma/r)^12 - (sigma/r)^6 ]
 *
 * @author Federico Celauro
 * @date   Apr 2, 2026
 */

#ifndef LJPART_H_
#define LJPART_H_

#include <cmath>
#include "../MassivePart/MassivePart.h"

/**
 * @class LJPart
 * @brief A MassivePart that stores Lennard-Jones force, potential, and
 *        interaction parameters (sigma, epsilon).
 *
 * This OpenMP variant adds the addForceAtomic() method for thread-safe
 * force accumulation using `#pragma omp atomic` directives.
 */
class LJPart : public MassivePart {
public:
    /// @name Constructors & Destructor
    /// @{
    LJPart();
    LJPart(long long int id, double m,
           double x, double y, double z,
           double Vx, double Vy, double Vz,
           double Ax, double Ay, double Az);
    LJPart(long long int id, double m,
           Eigen::Vector3d position, Eigen::Vector3d velocity,
           Eigen::Vector3d acceleration);
    LJPart(long long int id, double m,
           double x, double y, double z,
           double Vx, double Vy, double Vz);
    LJPart(long long int id, double m,
           Eigen::Vector3d position, Eigen::Vector3d velocity);
    LJPart(std::string dataString, char delimiter);
    virtual ~LJPart();
    /// @}

    /// @name Force component accessors (primary)
    /// @{
    double getForceX() const;
    void setForceX(double fx);
    double getForceY() const;
    void setForceY(double fy);
    double getForceZ() const;
    void setForceZ(double fz);
    /// @}

    /// @name Force component accessors (auxiliary)
    /// @{
    double getForceXAux() const;
    void setForceXAux(double fx);
    double getForceYAux() const;
    void setForceYAux(double fy);
    double getForceZAux() const;
    void setForceZAux(double fz);
    /// @}

    /// @name Force vector accessors
    /// @{
    const Eigen::Vector3d getForce() const;
    void setForce(const Eigen::Vector3d& force);
    const Eigen::Vector3d getForceAux() const;
    void setForceAux(const Eigen::Vector3d& force);
    /// @}

    /// @name LJ parameter accessors
    /// @{
    double getU() const;
    void setU(double u);
    double getSigma() const;
    void setSigma(double sigma);
    double getEpsilon() const;
    void setEpsilon(double epsilon);
    /// @}

    /// Sets the acceleration from the force
    void setAccelerationFromForce(void);

    /// Returns the magnitude of the primary force vector.
    double ForceMagnitude(void);

    /// Serializes the full particle state to a delimited string.
    std::string returnStateString(char delimiter) const;

    /// @name Pairwise LJ interaction computations
    /// @{
    Eigen::Vector3d computeForceij(const LJPart& other);
    double computeUij(const LJPart& other);
    Eigen::Vector3d computeForceijAndU(const LJPart& other, double& vij);
    /// @}

    /**
     * @brief Thread-safe force accumulation using OpenMP atomics.
     *
     * Used by the parallel force computation loop to safely update
     * the force vector from multiple threads simultaneously.
     */
    void addForceAtomic(const Eigen::Vector3d& force);

private:
    double U;                    ///< Per-particle potential energy
    Eigen::Vector3d f, fAux;     ///< Primary and auxiliary force vectors
    double sigma = 1.0, epsilon = 1.0;       ///< LJ zero-crossing distance and well depth

    /// Initializes force from acceleration; zeroes auxiliary force and potential.
    void initDefaultForce(void);

    /// NVI implementation for state serialization.
    virtual std::string returnStateStringPriv(char delimiter) const;

    /// @name Private force computation implementations
    /// @{
    Eigen::Vector3d computeForceijPriv(const LJPart& other);
    double computeUijPriv(const LJPart& other);
    Eigen::Vector3d computeForceijAndUPriv(const LJPart& other, double& vij);
    void addForceAtomicPriv(const Eigen::Vector3d& force);
    /// @}
};

#endif /* LJPART_H_ */
