/**
 * @file LagrangianPart.cpp
 * @brief Implementation of the LagrangianPart class.
 *
 * @author Federico Celauro
 * @date   Mar 4, 2026
 */

#include "../../include/physics/LagrangianPart/LagrangianPart.h"

// ---------------------------------------------------------------------------
// Constructors & Destructor
// ---------------------------------------------------------------------------

LagrangianPart::LagrangianPart() : PointInSpace() {
    initDefaultVelocity(0, 0, 0);
}

LagrangianPart::LagrangianPart(long long int id, double x, double y, double z,
                               double vx, double vy, double vz)
    : PointInSpace(id, x, y, z) {
    initDefaultVelocity(vx, vy, vz);
}

LagrangianPart::LagrangianPart(long long int id, Eigen::Vector3d position,
                               Eigen::Vector3d velocity)
    : PointInSpace(id, position) {
    initDefaultVelocity(velocity);
}

LagrangianPart::LagrangianPart(std::string allocator, char delimiter) {
    initFromString(allocator, delimiter);
}

LagrangianPart::~LagrangianPart() {
}

// ---------------------------------------------------------------------------
// Velocity component accessors
// ---------------------------------------------------------------------------

double LagrangianPart::getVx() const {
    return v.x();
}

void LagrangianPart::setVx(double vx) {
    v.x() = vx;
}

double LagrangianPart::getVy() const {
    return v.y();
}

void LagrangianPart::setVy(double vy) {
    v.y() = vy;
}

double LagrangianPart::getVz() const {
    return v.z();
}

void LagrangianPart::setVz(double vz) {
    v.z() = vz;
}

// ---------------------------------------------------------------------------
// Velocity vector accessors
// ---------------------------------------------------------------------------

const Eigen::Vector3d& LagrangianPart::getVelocity() const {
    return v;
}

void LagrangianPart::setVelocity(const Eigen::Vector3d& velocity) {
    v = velocity;
}

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------

double LagrangianPart::speed(void) {
    return v.norm();
}

// ---------------------------------------------------------------------------
// Initialization helpers
// ---------------------------------------------------------------------------

void LagrangianPart::initDefaultVelocity(double vx, double vy, double vz) {
    setVx(vx);
    setVy(vy);
    setVz(vz);
}

void LagrangianPart::initDefaultVelocity(Eigen::Vector3d& velocity) {
    setVelocity(velocity);
}

/**
 * @brief Parses a delimited string to initialize the particle.
 *
 * Expected format: "id, x, y, z, vx, vy, vz" (7 tokens).
 */
void LagrangianPart::initFromString(std::string dataString, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(dataString);

    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }

    if (tokens.size() != 7) {
        throw std::invalid_argument("Expected 7 CSV fields.");
    }

    setId(std::stoll(tokens[0]));
    setX(std::stod(tokens[1]));
    setY(std::stod(tokens[2]));
    setZ(std::stod(tokens[3]));

    initDefaultVelocity(
        std::stod(tokens[4]),
        std::stod(tokens[5]),
        std::stod(tokens[6])
    );
}

// ---------------------------------------------------------------------------
// State serialization (NVI pattern)
// ---------------------------------------------------------------------------

std::string LagrangianPart::returnStateString(char delimiter) const {
    return returnStateStringPriv(delimiter);
}

std::string LagrangianPart::returnStateStringPriv(char delimiter) const {
    std::stringstream ss;
    ss << getId()  << delimiter
       << getX()   << delimiter
       << getY()   << delimiter
       << getZ()   << delimiter
       << getVx()  << delimiter
       << getVy()  << delimiter
       << getVz();
    return ss.str();
}

// ---------------------------------------------------------------------------
// Position update
// ---------------------------------------------------------------------------

void LagrangianPart::updatePos(double dt) {
    updatePosPriv(dt);
}

/**
 * @brief Euler-step position update: r += v * dt.
 */
void LagrangianPart::updatePosPriv(double dt) {
    setPosition(getPosition() + v * dt);
}
