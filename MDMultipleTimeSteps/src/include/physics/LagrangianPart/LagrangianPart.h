/**
 * @file LagrangianPart.h
 * @brief Extends PointInSpace with a velocity vector for Lagrangian tracking.
 *
 * A LagrangianPart is a moving point: it stores position (inherited) and
 * velocity, and can update its position via x += v * dt.
 *
 * @author Federico Celauro
 * @date   Mar 4, 2026
 */

#ifndef LAGRANGIANPART_H_
#define LAGRANGIANPART_H_

#include "../../math/geometry/PointInSpace/PointInSpace.h"

/**
 * @class LagrangianPart
 * @brief A PointInSpace extended with a 3D velocity vector.
 *
 * Inherits the ID and position from PointInSpace and adds velocity
 * components.  Provides a simple Euler-step position update.
 */
class LagrangianPart : public PointInSpace {
public:
    /// @name Constructors & Destructor
    /// @{
    LagrangianPart();
    LagrangianPart(long long int id, double x, double y, double z,
                   double vx, double vy, double vz);
    LagrangianPart(long long int id, Eigen::Vector3d position,
                   Eigen::Vector3d velocity);
    LagrangianPart(std::string dataString, char delimiter);
    virtual ~LagrangianPart();
    /// @}

    /// @name Velocity component accessors
    /// @{
    double getVx() const;
    void setVx(double vx);

    double getVy() const;
    void setVy(double vy);

    double getVz() const;
    void setVz(double vz);
    /// @}

    /// @name Velocity vector accessors
    /// @{
    const Eigen::Vector3d& getVelocity() const;
    void setVelocity(const Eigen::Vector3d& velocity);
    /// @}

    /// Returns the magnitude of the velocity (speed).
    double speed(void);

    /// Serializes the state to a delimited string.
    std::string returnStateString(char delimiter) const;

    /// Advances the position by v * dt.
    void updatePos(double dt);

private:
    Eigen::Vector3d v;  ///< Velocity vector (vx, vy, vz)

    /// @name Initialization helpers
    /// @{
    void initDefaultVelocity(double vx, double vy, double vz);
    void initDefaultVelocity(Eigen::Vector3d& velocity);
    virtual void initFromString(std::string dataString, char delimiter);
    /// @}

    /// NVI implementation for state serialization.
    virtual std::string returnStateStringPriv(char delimiter) const;

    /// NVI implementation for position update.
    void updatePosPriv(double dt);
};

#endif /* LAGRANGIANPART_H_ */
