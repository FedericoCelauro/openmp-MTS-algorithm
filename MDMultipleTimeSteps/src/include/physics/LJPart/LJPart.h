/**
 * @file LJPart.h
 * @brief Extends MassivePart with Lennard-Jones interaction forces.
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
 * Maintains two force vectors:
 *  - @c f    : the total (or primary) force acting on the particle.
 *  - @c fAux : an auxiliary force buffer used in multiple-time-step schemes.
 *
 * Setting the force also automatically updates the acceleration
 * via Newton's second law (a = F / m).
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

    /// Returns the magnitude of the primary force vector.
    double ForceMagnitude(void);

    /// Serializes the full particle state to a delimited string.
    std::string returnStateString(char delimiter) const;

    /// @name Pairwise LJ interaction computations
    /// @{

    /// Computes the LJ force exerted by @p other on this particle.
    Eigen::Vector3d computeForceij(const LJPart& other);

    /// Computes the LJ pair potential U_ij with @p other.
    double computeUij(const LJPart& other);

    /**
     * @brief Computes both force and potential in a single pass.
     * @param other  The other particle.
     * @param[out] vij  Pair potential energy (overwritten).
     * @return The force vector on this particle due to @p other.
     */
    Eigen::Vector3d computeForceijAndU(const LJPart& other, double& vij);
    /// @}

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
    /// @}
};

#endif /* LJPART_H_ */
