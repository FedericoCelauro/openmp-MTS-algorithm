/**
 * @file PointInSpace.cpp
 * @brief Implementation of the PointInSpace class.
 *
 * @author Federico Celauro
 * @date   Mar 4, 2026
 */

#include "../../../include/math/geometry/PointInSpace/PointInSpace.h"

// ---------------------------------------------------------------------------
// Constructors & Destructor
// ---------------------------------------------------------------------------

PointInSpace::PointInSpace() {
    initDefaultPosition(0, 0, 0, 0);
}

PointInSpace::PointInSpace(long long int id, double x, double y, double z) {
    initDefaultPosition(id, x, y, z);
}

PointInSpace::PointInSpace(long long int id, Eigen::Vector3d position) {
    initDefaultPosition(id, position);
}

PointInSpace::PointInSpace(std::string allocator, char delimiter) {
    initFromString(allocator, delimiter);
}

PointInSpace::~PointInSpace() {
    // Eigen::Vector3d is stack-allocated; nothing to clean up.
}

// ---------------------------------------------------------------------------
// ID accessors
// ---------------------------------------------------------------------------

long long int PointInSpace::getId() const {
    return id;
}

void PointInSpace::setId(long long int id) {
    this->id = id;
}

// ---------------------------------------------------------------------------
// Coordinate accessors
// ---------------------------------------------------------------------------

double PointInSpace::getX() const {
    return r.x();
}

void PointInSpace::setX(double x) {
    this->r.x() = x;
}

double PointInSpace::getY() const {
    return r.y();
}

void PointInSpace::setY(double y) {
    this->r.y() = y;
}

double PointInSpace::getZ() const {
    return r.z();
}

void PointInSpace::setZ(double z) {
    this->r.z() = z;
}

// ---------------------------------------------------------------------------
// Position vector accessors
// ---------------------------------------------------------------------------

const Eigen::Vector3d& PointInSpace::getPosition() const {
    return r;
}

void PointInSpace::setPosition(const Eigen::Vector3d& position) {
    r = position;
}

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------

double PointInSpace::norm(void) const {
    return r.norm();
}

// ---------------------------------------------------------------------------
// Initialization helpers
// ---------------------------------------------------------------------------

void PointInSpace::initDefaultPosition(long long int id, double x, double y, double z) {
    setId(id);
    setX(x);
    setY(y);
    setZ(z);
}

void PointInSpace::initDefaultPosition(long long int id, Eigen::Vector3d& position) {
    setId(id);
    setPosition(position);
}

/**
 * @brief Parses a delimited string to initialise the point.
 *
 * Expected format: "id, x, y, z" (4 tokens separated by @p delimiter).
 */
void PointInSpace::initFromString(std::string dataString, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(dataString);

    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }

    if (tokens.size() != 4) {
        throw std::invalid_argument("Expected 4 CSV fields.");
    }

    initDefaultPosition(
        std::stoll(tokens[0]),
        std::stod(tokens[1]),
        std::stod(tokens[2]),
        std::stod(tokens[3])
    );
}

// ---------------------------------------------------------------------------
// State serialization (NVI pattern)
// ---------------------------------------------------------------------------

std::string PointInSpace::returnStateString(char delimiter) const {
    return returnStateStringPriv(delimiter);
}

std::string PointInSpace::returnStateStringPriv(char delimiter) const {
    std::stringstream ss;
    ss << getId() << delimiter
       << getX()  << delimiter
       << getY()  << delimiter
       << getZ();
    return ss.str();
}
