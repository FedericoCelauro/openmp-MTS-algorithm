/**
 * @file LJPart.cpp
 * @brief Implementation of the LJPart class (Lennard-Jones particle).
 *
 * All pairwise interaction formulas follow the standard 12-6 LJ potential:
 *   U(r) = 4 * epsilon * [ (sigma/r)^12 - (sigma/r)^6 ]
 *
 * The force is derived analytically:
 *   F(r) = 24 * epsilon * [ 2*(sigma/r)^12 - (sigma/r)^6 ] / r^2  *  r_vec
 *
 * @author Federico Celauro
 * @date   Apr 2, 2026
 */

#include "../../include/physics/LJPart/LJPart.h"
#include <iomanip>
#include <limits>

// ---------------------------------------------------------------------------
// Constructors & Destructor
// ---------------------------------------------------------------------------

LJPart::LJPart() : MassivePart() {
    initDefaultForce();
}

LJPart::LJPart(long long int id, double m,
               double x, double y, double z,
               double Vx, double Vy, double Vz,
               double Ax, double Ay, double Az)
    : MassivePart(id, m, x, y, z, Vx, Vy, Vz, Ax, Ay, Az) {
    initDefaultForce();
}

LJPart::LJPart(long long int id, double m,
               Eigen::Vector3d position, Eigen::Vector3d velocity,
               Eigen::Vector3d acceleration)
    : MassivePart(id, m, position, velocity, acceleration) {
    initDefaultForce();
}

LJPart::LJPart(long long int id, double m,
               double x, double y, double z,
               double Vx, double Vy, double Vz)
    : MassivePart(id, m, x, y, z, Vx, Vy, Vz) {
    initDefaultForce();
}

LJPart::LJPart(long long int id, double m,
               Eigen::Vector3d position, Eigen::Vector3d velocity)
    : MassivePart(id, m, position, velocity) {
    initDefaultForce();
}

LJPart::LJPart(std::string allocator, char delimiter)
    : MassivePart(allocator, delimiter) {
    initDefaultForce();
}

LJPart::~LJPart() {
}

// ---------------------------------------------------------------------------
// Force component accessors (primary)
// ---------------------------------------------------------------------------

double LJPart::getForceX() const {
    return f.x();
}

void LJPart::setForceX(double forceX) {
    this->f.x() = forceX;
    this->setAx(forceX / this->getMass());
}

double LJPart::getForceY() const {
    return f.y();
}

void LJPart::setForceY(double forceY) {
    this->f.y() = forceY;
    this->setAy(forceY / this->getMass());
}

double LJPart::getForceZ() const {
    return f.z();
}

void LJPart::setForceZ(double forceZ) {
    this->f.z() = forceZ;
    this->setAz(forceZ / this->getMass());
}

// ---------------------------------------------------------------------------
// Force component accessors (auxiliary)
// ---------------------------------------------------------------------------

double LJPart::getForceXAux() const {
    return fAux.x();
}

void LJPart::setForceXAux(double forceX) {
    this->fAux.x() = forceX;
}

double LJPart::getForceYAux() const {
    return fAux.y();
}

void LJPart::setForceYAux(double forceY) {
    this->fAux.y() = forceY;
}

double LJPart::getForceZAux() const {
    return fAux.z();
}

void LJPart::setForceZAux(double forceZ) {
    this->fAux.z() = forceZ;
}

// ---------------------------------------------------------------------------
// Force vector accessors
// ---------------------------------------------------------------------------

const Eigen::Vector3d LJPart::getForce() const {
    return f;
}

void LJPart::setForce(const Eigen::Vector3d& force) {
    f = force;
    this->setAcceleration(force / this->getMass());
}

const Eigen::Vector3d LJPart::getForceAux() const {
    return fAux;
}

void LJPart::setForceAux(const Eigen::Vector3d& force) {
    fAux = force;
}

// ---------------------------------------------------------------------------
// LJ parameter accessors
// ---------------------------------------------------------------------------

double LJPart::getU() const {
    return U;
}

void LJPart::setU(double u) {
    this->U = u;
}

double LJPart::getSigma() const {
    return sigma;
}

void LJPart::setSigma(double sigma) {
    this->sigma = sigma;
}

double LJPart::getEpsilon() const {
    return epsilon;
}

void LJPart::setEpsilon(double epsilon) {
    this->epsilon = epsilon;
}

// ---------------------------------------------------------------------------
// Initialization
// ---------------------------------------------------------------------------

void LJPart::initDefaultForce(void) {
    setForce(getAcceleration() * getMass());
    setForceXAux(0);
    setForceYAux(0);
    setForceZAux(0);
    setU(0);
}

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------

double LJPart::ForceMagnitude(void) {
    return f.norm();
}

// ---------------------------------------------------------------------------
// State serialization (NVI pattern)
// ---------------------------------------------------------------------------

std::string LJPart::returnStateString(char delimiter) const {
    return returnStateStringPriv(delimiter);
}

std::string LJPart::returnStateStringPriv(char delimiter) const {
    std::stringstream ss;
    ss << std::setprecision(std::numeric_limits<double>::max_digits10)
       << getId()         << delimiter
       << getMass()       << delimiter
       << getX()          << delimiter
       << getY()          << delimiter
       << getZ()          << delimiter
       << getVx()         << delimiter
       << getVy()         << delimiter
       << getVz()         << delimiter
       << getAx()         << delimiter
       << getAy()         << delimiter
       << getAz()         << delimiter
       << getForceX()     << delimiter
       << getForceY()     << delimiter
       << getForceZ()     << delimiter
       << getForceXAux()  << delimiter
       << getForceYAux()  << delimiter
       << getForceZAux()  << delimiter
       << getU();
    return ss.str();
}

// ---------------------------------------------------------------------------
// Pairwise LJ interaction computations
// ---------------------------------------------------------------------------

Eigen::Vector3d LJPart::computeForceij(const LJPart& particle) {
    return computeForceijPriv(particle);
}

double LJPart::computeUij(const LJPart& particle) {
    return computeUijPriv(particle);
}

Eigen::Vector3d LJPart::computeForceijAndU(const LJPart& particle, double& vij) {
    return computeForceijAndUPriv(particle, vij);
}

/**
 * @brief Computes the LJ force between this particle and @p particle.
 *
 * Uses the standard 12-6 LJ force formula with no cutoff applied.
 * F_ij = 24*eps * [2*(sigma/r)^12 - (sigma/r)^6] / r^2 * r_vec
 */
Eigen::Vector3d LJPart::computeForceijPriv(const LJPart& particle) {
    Eigen::Vector3d rij = getPosition() - particle.getPosition();
    const double distance = rij.norm();

    double multFactor = 24 * epsilon
        * (2 * std::pow(sigma / distance, 12)
              - std::pow(sigma / distance, 6))
        / (distance * distance);

    return rij * multFactor;
}

/**
 * @brief Computes the LJ pair potential between this particle and @p particle.
 *
 * U_ij = 4*eps * [(sigma/r)^12 - (sigma/r)^6]
 */
double LJPart::computeUijPriv(const LJPart& particle) {
    double distance = (getPosition() - particle.getPosition()).norm();
    return 4.0 * epsilon
        * (std::pow(sigma / distance, 12) - std::pow(sigma / distance, 6));
}

/**
 * @brief Computes both the LJ force and potential in a single pass.
 *
 * Avoids redundant distance calculations by using the squared distance
 * and reduced variables: srij2 = (sigma/r)^2, srij6 = srij2^3, etc.
 *
 * @param particle  The other particle.
 * @param[out] vij  The pair potential energy (overwritten).
 * @return The force vector on this particle due to @p particle.
 */
Eigen::Vector3d LJPart::computeForceijAndUPriv(const LJPart& particle, double& vij) {
    Eigen::Vector3d rij = getPosition() - particle.getPosition();

    const double rij2   = rij.squaredNorm();
    const double srij2  = sigma * sigma / rij2;
    const double srij6  = srij2 * srij2 * srij2;
    const double srij12 = srij6 * srij6;

    vij = 4 * epsilon * (srij12 - srij6);
    double multFactor = 24 * epsilon * (2 * srij12 - srij6) / rij2;
    return rij * multFactor;
}
