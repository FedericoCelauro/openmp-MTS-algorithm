/**
 * @file PointInSpace.h
 * @brief Defines a point in 3D Euclidean space with an identifier.
 *
 * This is the base class for all particle types in the simulation.
 * It stores a unique particle ID and a 3D position vector using Eigen.
 *
 * @author Federico Celauro
 * @date   Mar 4, 2026
 */

#ifndef POINTINSPACE_H_
#define POINTINSPACE_H_

#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

/**
 * @class PointInSpace
 * @brief Represents a point in 3D space identified by a unique ID.
 *
 * Provides getters/setters for individual coordinates (x, y, z) and the
 * full Eigen::Vector3d position vector.  Supports initialization from
 * constructor arguments or from a delimited string (CSV parsing).
 */
class PointInSpace {

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /// @name Constructors & Destructor
    /// @{
    PointInSpace();
    PointInSpace(long long int id, double x, double y, double z);
    PointInSpace(long long int id, Eigen::Vector3d position);
    PointInSpace(std::string dataString, char delimiter);
    virtual ~PointInSpace();
    /// @}

    /// @name ID accessors
    /// @{
    long long int getId() const;
    void setId(long long int id);
    /// @}

    /// @name Coordinate accessors
    /// @{
    double getX() const;
    void setX(double x);

    double getY() const;
    void setY(double y);

    double getZ() const;
    void setZ(double z);
    /// @}

    /// @name Position vector accessors
    /// @{
    const Eigen::Vector3d& getPosition() const;
    void setPosition(const Eigen::Vector3d& position);
    /// @}

    /// Returns the Euclidean distance from the origin.
    double norm(void) const;

    /// Serializes the point state to a delimited string.
    std::string returnStateString(char delimiter) const;

private:
    long long int id;   ///< Unique particle identifier
    Eigen::Vector3d r;  ///< Position vector (x, y, z)

    /// @name Initialization helpers
    /// @{
    void initDefaultPosition(long long int id, double x, double y, double z);
    void initDefaultPosition(long long int id, Eigen::Vector3d& position);
    virtual void initFromString(std::string dataString, char delimiter);
    /// @}

    /// Formats the internal state as a delimited string (NVI pattern).
    virtual std::string returnStateStringPriv(char delimiter) const;
};

#endif /* POINTINSPACE_H_ */
