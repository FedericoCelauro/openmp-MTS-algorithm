/**
 * @file SimParams.h
 * @brief Aggregates all simulation parameters for a Molecular Dynamics run.
 *
 * Reads physical and numerical parameters from a configuration file and
 * derives box dimensions, particle counts, and other simulation constants.
 *
 * @author Federico Celauro
 * @date   Apr 6, 2026
 */

#ifndef INCLUDE_SIMPARAMS_SIMPARAMS_H_
#define INCLUDE_SIMPARAMS_SIMPARAMS_H_

#include <string>
#include <map>
#include <any>

#include "../IO/ASCII/asciiFileReader.h"

/**
 * @class SimParams
 * @brief Container for all simulation parameters (physical, numerical,
 *        and derived quantities).
 *
 * Can be constructed with default values or by parsing a key-value
 * configuration file. Derived quantities are computed during construction.
 * After calling setters, recompute the corresponding derived quantities explicitly.
 */
class SimParams {
public:
    SimParams();

    /**
     * @brief Constructs SimParams by reading a configuration file.
     * @param path        Directory containing the config file.
     * @param configFile  Name of the configuration file.
     */
    SimParams(std::string path, std::string configFile);

    virtual ~SimParams();

    /// @name Time-stepping parameters
    /// @{
    double getDt() const;
    void setDt(double dt);
    long int getTotsteps() const;
    void setTotsteps(long int totsteps);
    double getTotTime() const;
    void setTotTime(void);   ///< Computes totTime = totsteps * dt
    /// @}

    /// @name LJ interaction parameters
    /// @{
    double getEpsilon() const;
    void setEpsilon(double epsilon);
    double getSigma() const;
    void setSigma(double sigma);
    double getRc() const;
    void setRc(double rc);
    /// @}

    /// @name Particle masses
    /// @{
    double getMass1() const;
    void setMass1(double mass1);
    double getMass2() const;
    void setMass2(double mass2);
    double getNPartRatio() const;
    void setNPartRatio(double nratio);
    /// @}

    /// @name Thermodynamic parameters
    /// @{
    double getT() const;
    void setT(double t);
    double getRho() const;
    void setRho(double rho);
    /// @}

    /// @name Grid & derived geometry
    /// @{
    int getNx() const;
    void setNx(int nx);
    int getNy() const;
    void setNy(int ny);
    int getNz() const;
    void setNz(int nz);
    double getLx() const;
    void setLx(void);   ///< Derives Lx from nx, sigma, and rho
    double getLy() const;
    void setLy(void);   ///< Derives Ly from ny, sigma, and rho
    double getLz() const;
    void setLz(void);   ///< Derives Lz from nz, sigma, and rho
    /// @}

    /// @name Particle count accessors
    /// @{
    long int getNparticles() const;
    void setNparticles(void);    ///< Computes Nparticles = nx * ny * nz
    long int getNPartSmallMass() const;
    void setNPartSmallMass(void);
    int getRatio() const;
    void setRatio(void);         ///< Computes mass ratio sqrt(m_heavy / m_light)
    /// @}

    /// @name Output path accessors
    /// @{
    const std::string& getOutputTimePath() const;
    void setOutputTimePath(const std::string& outputTimePath);
    const std::string& getEnergyPath() const;
    void setEnergyPath(const std::string& energyPath);
    /// @}

    /// Prints all parameter values to stdout.
    void printSimParams(void);

private:
    // Grid dimensions
    int nx, ny, nz, Ratio;

    // Physical and numerical parameters
    double rc, sigma, epsilon, mass1, mass2, NPartRatio;
    double dt, T, rho, totTime, Lx, Ly, Lz;

    // Derived particle counts
    long int totsteps, Nparticles, NPartSmallMass;

    // Output paths
    std::string OutputTimePath, EnergyPath;

    /// Private implementation of printSimParams.
    void printSimParamsPriv(void);
};

#endif /* INCLUDE_SIMPARAMS_SIMPARAMS_H_ */
