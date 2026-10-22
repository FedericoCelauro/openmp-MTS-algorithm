/**
 * @file MDPart.h
 * @brief Extends LJPart with periodic boundary conditions and cutoff radius.
 *
 * An MDPart is an LJ particle placed inside a periodic simulation box
 * of dimensions L = (Lx, Ly, Lz).  It applies the minimum image
 * convention when computing interparticle distances and enforces
 * periodic boundary conditions (PBC) on positions.
 *
 * @author Federico Celauro
 * @date   Apr 6, 2026
 */

#ifndef MDPART_MDPART_H_
#define MDPART_MDPART_H_

#include "../physics/LJPart/LJPart.h"

/**
 * @class MDPart
 * @brief LJPart with a periodic simulation box and interaction cutoff.
 *
 * Forces and potentials are computed only for pairs within the cutoff
 * distance @c rc.  A constant potential shift (Ucut) is subtracted from the
 * potential to ensure continuity at the cutoff.
 */
class MDPart : public LJPart {
public:
    /// @name Constructors & Destructor
    /// @{
    MDPart();
    MDPart(long long int id, double m,
           double x, double y, double z,
           double Vx, double Vy, double Vz,
           double Ax, double Ay, double Az,
           double Lx, double Ly, double Lz, double rc);
    MDPart(long long int id, double m,
           Eigen::Vector3d position, Eigen::Vector3d velocity,
           Eigen::Vector3d acceleration, Eigen::Vector3d l, double rc);
    MDPart(long long int id, double m,
           double x, double y, double z,
           double Vx, double Vy, double Vz,
           double Ax, double Ay, double Az);
    MDPart(long long int id, double m,
           Eigen::Vector3d position, Eigen::Vector3d velocity,
           Eigen::Vector3d acceleration, double rc);
    MDPart(long long int id, double m,
           double x, double y, double z,
           double Vx, double Vy, double Vz);
    MDPart(long long int id, double m,
           Eigen::Vector3d position, Eigen::Vector3d velocity);
    MDPart(std::string dataString, char delimiter);
    virtual ~MDPart();
    /// @}

    /// @name Box dimension accessors
    /// @{
    double getLx() const;
    void setLx(double lx);
    double getLy() const;
    void setLy(double ly);
    double getLz() const;
    void setLz(double lz);
    Eigen::Vector3d getL() const;
    void setL(Eigen::Vector3d l);
    /// @}

    /// @name Cutoff accessors
    /// @{
    double getRc() const;
    void setRc(double rc);
    double getUcut() const;
    void setUcut(void);  ///< Computes Ucut from current sigma, epsilon, rc
    /// @}

    /// Applies periodic boundary conditions (wraps position into the box).
    void pbc();

    /// @name Pairwise interaction computations (with cutoff & minimum image)
    /// @{
    Eigen::Vector3d computeForceij(const MDPart& other);
    double computeUij(const MDPart& other);

    /**
     * @brief Computes both force and potential in a single pass.
     * @param other  The other particle.
     * @param[in,out] vij  Pair potential energy (accumulated, not overwritten).
     * @return Force vector on this particle due to @p other.
     */
    Eigen::Vector3d computeForceijAndU(const MDPart& other, double& vij);
    /// @}

private:
    Eigen::Vector3d L;   ///< Simulation box dimensions (Lx, Ly, Lz)
    double rc;           ///< Interaction cutoff radius
    double Ucut = 0.0;         ///< LJ potential value at the cutoff (for shifting)

    /// @name Private force computation implementations
    /// @{
    Eigen::Vector3d computeForceijPriv(const MDPart& other);
    double computeUijPriv(const MDPart& other);
    Eigen::Vector3d computeForceijAndUPriv(const MDPart& other, double& vij);
    /// @}
};

#endif /* MDPART_MDPART_H_ */
