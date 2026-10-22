/**
 * @file MassivePart.h
 * @brief Extends LagrangianPart with mass and acceleration.
 *
 * A MassivePart is a Lagrangian particle that also carries a mass and
 * an acceleration vector.  It can update its velocity via v += a * dt.
 *
 * @author Federico Celauro
 * @date   Mar 4, 2026
 */

#ifndef MASSIVEPART_H_
#define MASSIVEPART_H_

#include "../LagrangianPart/LagrangianPart.h"
#include <stdexcept>

/**
 * @class MassivePart
 * @brief A LagrangianPart extended with mass and 3D acceleration.
 *
 * Inherits position and velocity from LagrangianPart and adds mass
 * (with a positivity check) and acceleration components.
 */
class MassivePart : public LagrangianPart {
public:
    /// @name Constructors & Destructor
    /// @{
    MassivePart();
    MassivePart(long long int id, double m,
                double x, double y, double z,
                double Vx, double Vy, double Vz,
                double Ax, double Ay, double Az);
    MassivePart(long long int id, double m,
                Eigen::Vector3d position, Eigen::Vector3d velocity,
                Eigen::Vector3d acceleration);
    MassivePart(long long int id, double m,
                double x, double y, double z,
                double Vx, double Vy, double Vz);
    MassivePart(long long int id, double m,
                Eigen::Vector3d position, Eigen::Vector3d velocity);
    MassivePart(std::string dataString, char delimiter);
    virtual ~MassivePart();
    /// @}

    /// @name Mass accessors
    /// @{
    double getMass() const;
    void setMass(double mass);  ///< @throws std::invalid_argument if mass is non-finite or <= 0
    /// @}

    /// @name Acceleration component accessors
    /// @{
    double getAx() const;
    void setAx(double ax);

    double getAy() const;
    void setAy(double ay);

    double getAz() const;
    void setAz(double az);
    /// @}

    /// @name Acceleration vector accessors
    /// @{
    const Eigen::Vector3d& getAcceleration() const;
    void setAcceleration(const Eigen::Vector3d& acceleration);
    /// @}

    /// Returns the magnitude of the acceleration vector.
    double accMagnitude(void);

    /// Serializes the state to a delimited string.
    virtual std::string returnStateString(char delimiter) const;

    /// Advances the velocity by a * dt.
    void updateVel(double dt);

    /**
     * @brief Thread-safe force accumulation using OpenMP atomics.
     *
     * Used by the parallel force computation loop to safely update
     * the acceleration vector from multiple threads simultaneously.
     */
    void addAccelerationAtomic(const Eigen::Vector3d&);


private:
    double mass;        ///< Particle mass (must be positive)
    Eigen::Vector3d a = Eigen::Vector3d::Zero();  ///< Acceleration vector (ax, ay, az)

    /// @name Initialization helpers
    /// @{
    void initDefaultMass(double m);
    void initDefaultAcceleration(double Ax, double Ay, double Az);
    void initDefaultAcceleration(Eigen::Vector3d& acceleration);
    virtual void initFromString(std::string dataString, char delimiter);
    /// @}

    /// NVI implementation for state serialization.
    virtual std::string returnStateStringPriv(char delimiter) const;

    /// NVI implementation for velocity update.
    void updateVelPriv(double dt);

    /// NVI implementation for multiple-thread-safe acceleration update.
    void addAccelerationAtomicPriv(const Eigen::Vector3d&);

};

#endif /* MASSIVEPART_H_ */
