/**
 * @file MassivePart.cpp
 * @brief Implementation of the MassivePart class.
 *
 * @author Federico Celauro
 * @date   Mar 4, 2026
 */

#include "../../include/physics/MassivePart/MassivePart.h"
#include <cmath>

// ---------------------------------------------------------------------------
// Constructors & Destructor
// ---------------------------------------------------------------------------

MassivePart::MassivePart() : LagrangianPart() {
    initDefaultAcceleration(0, 0, 0);
    initDefaultMass(1);
}

MassivePart::MassivePart(long long int id, double m,
                         double x, double y, double z,
                         double Vx, double Vy, double Vz,
                         double Ax, double Ay, double Az)
    : LagrangianPart(id, x, y, z, Vx, Vy, Vz) {
    initDefaultAcceleration(Ax, Ay, Az);
    initDefaultMass(m);
}

MassivePart::MassivePart(long long int id, double m,
                         Eigen::Vector3d position,
                         Eigen::Vector3d velocity,
                         Eigen::Vector3d acceleration)
    : LagrangianPart(id, position, velocity) {
    initDefaultAcceleration(acceleration);
    initDefaultMass(m);
}

MassivePart::MassivePart(long long int id, double m,
                         double x, double y, double z,
                         double Vx, double Vy, double Vz)
    : LagrangianPart(id, x, y, z, Vx, Vy, Vz) {
    initDefaultMass(m);
}

MassivePart::MassivePart(long long int id, double m,
                         Eigen::Vector3d position,
                         Eigen::Vector3d velocity)
    : LagrangianPart(id, position, velocity) {
    initDefaultMass(m);
}

MassivePart::MassivePart(std::string allocator, char delimiter) {
    initFromString(allocator, delimiter);
}

MassivePart::~MassivePart() {
}

// ---------------------------------------------------------------------------
// Mass accessors
// ---------------------------------------------------------------------------

double MassivePart::getMass() const {
    return mass;
}

void MassivePart::setMass(double mass) {
    if (!std::isfinite(mass) || mass <= 0) {
        throw std::invalid_argument("Mass must be positive.");
    }
    this->mass = mass;
}

// ---------------------------------------------------------------------------
// Acceleration component accessors
// ---------------------------------------------------------------------------

double MassivePart::getAx() const {
    return a.x();
}

void MassivePart::setAx(double ax) {
    this->a.x() = ax;
}

double MassivePart::getAy() const {
    return a.y();
}

void MassivePart::setAy(double ay) {
    this->a.y() = ay;
}

double MassivePart::getAz() const {
    return a.z();
}

void MassivePart::setAz(double az) {
    this->a.z() = az;
}

// ---------------------------------------------------------------------------
// Acceleration vector accessors
// ---------------------------------------------------------------------------

const Eigen::Vector3d& MassivePart::getAcceleration() const {
    return a;
}

void MassivePart::setAcceleration(const Eigen::Vector3d& acceleration) {
    a = acceleration;
}

double MassivePart::accMagnitude(void) {
    return a.norm();
}

// ---------------------------------------------------------------------------
// Initialization helpers
// ---------------------------------------------------------------------------

void MassivePart::initDefaultMass(double m) {
    setMass(m);
}

void MassivePart::initDefaultAcceleration(double Ax, double Ay, double Az) {
    setAx(Ax);
    setAy(Ay);
    setAz(Az);
}

void MassivePart::initDefaultAcceleration(Eigen::Vector3d& acceleration) {
    setAcceleration(acceleration);
}

/**
 * @brief Parses a delimited string to initialize the particle.
 *
 * Expected format: "id, m, x, y, z, vx, vy, vz" (8 tokens).
 */
void MassivePart::initFromString(std::string dataString, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(dataString);

    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }

    if (tokens.size() != 8) {
        throw std::invalid_argument("Expected 8 CSV fields.");
    }

    for (size_t i = 0; i < tokens.size(); ++i) {
        size_t used = 0;
        if (i == 0) {
            std::stoll(tokens[i], &used);
        } else if (!std::isfinite(std::stod(tokens[i], &used))) {
            throw std::invalid_argument("CSV fields must be finite numbers.");
        }
        if (tokens[i].find_first_not_of(" \t\r", used) != std::string::npos) {
            throw std::invalid_argument("Invalid numeric CSV field: " + tokens[i]);
        }
    }
    setId(std::stoll(tokens[0]));
    initDefaultMass(std::stod(tokens[1]));
    setX(std::stod(tokens[2]));
    setY(std::stod(tokens[3]));
    setZ(std::stod(tokens[4]));
    setVx(std::stod(tokens[5]));
    setVy(std::stod(tokens[6]));
    setVz(std::stod(tokens[7]));
}

// ---------------------------------------------------------------------------
// State serialization (NVI pattern)
// ---------------------------------------------------------------------------

std::string MassivePart::returnStateString(char delimiter) const {
    return returnStateStringPriv(delimiter);
}

std::string MassivePart::returnStateStringPriv(char delimiter) const {
    std::stringstream ss;
    ss << getId()    << delimiter
       << getMass()  << delimiter
       << getX()     << delimiter
       << getY()     << delimiter
       << getZ()     << delimiter
       << getVx()    << delimiter
       << getVy()    << delimiter
       << getVz();
    return ss.str();
}

// ---------------------------------------------------------------------------
// Velocity update
// ---------------------------------------------------------------------------

void MassivePart::updateVel(double dt) {
    updateVelPriv(dt);
}

/**
 * @brief Euler-step velocity update: v += a * dt.
 */
void MassivePart::updateVelPriv(double dt) {
    setVelocity(getVelocity() + a * dt);
}
