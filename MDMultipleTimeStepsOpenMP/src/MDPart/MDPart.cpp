/**
 * @file MDPart.cpp
 * @brief Implementation of the MDPart class (OpenMP variant).
 *
 * All pairwise interaction computations apply the minimum image
 * convention before evaluating forces/potentials, and enforce the
 * cutoff radius rc.  The potential is shifted so that U(rc) = 0.
 *
 * @author Federico Celauro
 * @date   Apr 6, 2026
 */

#include "../include/MDPart/MDPart.h"

// ---------------------------------------------------------------------------
// Constructors & Destructor
// ---------------------------------------------------------------------------

MDPart::MDPart() : LJPart() {
    setLx(0); setLy(0); setLz(0);
    setRc(0);
}

MDPart::MDPart(long long int id, double m,
               double x, double y, double z,
               double Vx, double Vy, double Vz,
               double Ax, double Ay, double Az,
               double Lx, double Ly, double Lz, double rc)
    : LJPart(id, m, x, y, z, Vx, Vy, Vz, Ax, Ay, Az) {
    setLx(Lx); setLy(Ly); setLz(Lz);
    setRc(rc);
}

MDPart::MDPart(long long int id, double m,
               Eigen::Vector3d position, Eigen::Vector3d velocity,
               Eigen::Vector3d acceleration, Eigen::Vector3d l, double rc)
    : LJPart(id, m, position, velocity, acceleration) {
    setL(l);
    setRc(rc);
}

MDPart::MDPart(long long int id, double m,
               double x, double y, double z,
               double Vx, double Vy, double Vz,
               double Ax, double Ay, double Az)
    : LJPart(id, m, x, y, z, Vx, Vy, Vz, Ax, Ay, Az) {
    setLx(0); setLy(0); setLz(0);
    setRc(0);
}

MDPart::MDPart(long long int id, double m,
               Eigen::Vector3d position, Eigen::Vector3d velocity,
               Eigen::Vector3d acceleration, double rc)
    : LJPart(id, m, position, velocity, acceleration) {
    setLx(0); setLy(0); setLz(0);
    setRc(rc);
}

MDPart::MDPart(long long int id, double m,
               double x, double y, double z,
               double Vx, double Vy, double Vz)
    : LJPart(id, m, x, y, z, Vx, Vy, Vz) {
    setLx(0); setLy(0); setLz(0);
    setRc(0);
}

MDPart::MDPart(long long int id, double m,
               Eigen::Vector3d position, Eigen::Vector3d velocity)
    : LJPart(id, m, position, velocity) {
    setLx(0); setLy(0); setLz(0);
    setRc(0);
}

MDPart::MDPart(std::string allocator, char delimiter)
    : LJPart(allocator, delimiter) {
    setLx(0); setLy(0); setLz(0);
    setRc(0);
}

MDPart::~MDPart() {
}

// ---------------------------------------------------------------------------
// Box dimension accessors
// ---------------------------------------------------------------------------

double MDPart::getLx() const { return L.x(); }
void MDPart::setLx(double L) { this->L.x() = L; }

double MDPart::getLy() const { return L.y(); }
void MDPart::setLy(double L) { this->L.y() = L; }

double MDPart::getLz() const { return L.z(); }
void MDPart::setLz(double L) { this->L.z() = L; }

Eigen::Vector3d MDPart::getL() const { return L; }
void MDPart::setL(Eigen::Vector3d l) { L = l; }

// ---------------------------------------------------------------------------
// Cutoff accessors
// ---------------------------------------------------------------------------

double MDPart::getRc() const { return rc; }
void MDPart::setRc(double rc) { this->rc = rc; }

double MDPart::getUcut() const { return Ucut; }

/**
 * @brief Computes the LJ potential at the cutoff radius.
 *
 * Subtracting this constant makes U(rc) = 0 and the truncated potential
 * continuous. The force is not shifted and remains discontinuous at rc.
 */
void MDPart::setUcut() {
    const double rijsq  = rc * rc;
    const double srij2  = getSigma() * getSigma() / rijsq;
    const double srij6  = srij2 * srij2 * srij2;
    const double srij12 = srij6 * srij6;
    this->Ucut = 4.0 * getEpsilon() * (srij12 - srij6);
}

// ---------------------------------------------------------------------------
// Periodic boundary conditions
// ---------------------------------------------------------------------------

/**
 * @brief Wraps the particle position into the primary simulation cell.
 *
 * Uses the nearest-integer rounding method:
 *   r_i = r_i - L_i * round(r_i / L_i)
 */
void MDPart::pbc(void) {
    setPosition(getPosition().array()
                - L.array() * (getPosition().array() / L.array()).round());
}

// ---------------------------------------------------------------------------
// Pairwise interaction computations (with cutoff & minimum image)
// ---------------------------------------------------------------------------

Eigen::Vector3d MDPart::computeForceij(const MDPart& particle) {
    return computeForceijPriv(particle);
}

double MDPart::computeUij(const MDPart& particle) {
    return computeUijPriv(particle);
}

Eigen::Vector3d MDPart::computeForceijAndU(const MDPart& particle, double& vij) {
    return computeForceijAndUPriv(particle, vij);
}

/**
 * @brief Computes the LJ force with minimum image convention and cutoff.
 *
 * The separation vector is adjusted using the minimum image convention
 * before evaluating the force.  If the distance exceeds rc, zero force
 * is returned.
 */
Eigen::Vector3d MDPart::computeForceijPriv(const MDPart& particle) {
    Eigen::Vector3d rij = getPosition() - particle.getPosition();

    // Minimum image convention
    rij.array() -= L.array() * (rij.array() / L.array()).round();

    const double distance = rij.norm();

    double multFactor = 0;
    if (distance <= rc) {
        multFactor = 24 * getEpsilon()
            * (2 * std::pow(getSigma() / distance, 12)
                  - std::pow(getSigma() / distance, 6))
            / (distance * distance);
    }

    return rij * multFactor;
}

/**
 * @brief Computes the shifted LJ potential with minimum image convention.
 *
 * Returns U(r) - Ucut if r <= rc, or 0 otherwise.
 */
double MDPart::computeUijPriv(const MDPart& particle) {
    Eigen::Vector3d rij = getPosition() - particle.getPosition();

    // Minimum image convention
    rij.array() -= L.array() * (rij.array() / L.array()).round();

    const double distance = rij.norm();

    if (distance <= rc) {
        double v = 4.0 * getEpsilon()
            * (std::pow(getSigma() / distance, 12)
               - std::pow(getSigma() / distance, 6));
        return v - Ucut;
    }
    return 0.0;
}

/**
 * @brief Computes both force and shifted potential in a single pass.
 *
 * Uses squared-distance arithmetic and reduced variables for efficiency.
 * The pair potential is accumulated into @p vij (not overwritten).
 */
Eigen::Vector3d MDPart::computeForceijAndUPriv(const MDPart& particle, double& vij) {
    Eigen::Vector3d rij = getPosition() - particle.getPosition();

    // Minimum image convention
    rij.array() -= L.array() * (rij.array() / L.array()).round();

    const double rij2 = rij.squaredNorm();

    double multFactor = 0;
    if (rij2 <= rc * rc) {
        const double srij2  = getSigma() * getSigma() / rij2;
        const double srij6  = srij2 * srij2 * srij2;
        const double srij12 = srij6 * srij6;

        double v = 4 * getEpsilon() * (srij12 - srij6);
        multFactor = 24 * getEpsilon() * (2 * srij12 - srij6) / rij2;
        vij += v - Ucut;
    }
    return rij * multFactor;
}
