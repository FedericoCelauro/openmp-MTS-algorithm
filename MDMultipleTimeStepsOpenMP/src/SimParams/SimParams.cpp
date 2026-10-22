/**
 * @file SimParams.cpp
 * @brief Implementation of the SimParams class.
 *
 * Provides configuration parsing using a variant-based dispatch
 * (std::variant + std::visit) to assign typed parameter values from
 * a text configuration file.
 *
 * @author Federico Celauro
 * @date   Apr 6, 2026
 */

#include "../include/SimParams/SimParams.h"

#include <variant>
#include <sstream>
#include <cmath>
#include <set>
#include <algorithm>
#include <limits>

// ---------------------------------------------------------------------------
// Default constructor — sets reasonable LJ defaults
// ---------------------------------------------------------------------------

SimParams::SimParams() {
    setNx(4); setNy(4); setNz(4);
    setRho(0.5);
    setRc(2.5);
    setSigma(1.0);
    setEpsilon(1.0);
    setTotsteps(800);
    setMass1(1.0);
    setMass2(100.0);
    setNPartRatio(0.01);
    setDt(0.05);
    setT(1.2);

    // Derive dependent quantities
    setNparticles();
    setNPartSmallMass();
    setTotTime();
    setRatio();
    setLx();
    setLy();
    setLz();

    setOutputTimePath("../jobs/sim1/OutputTimes");
    setEnergyPath("../jobs/sim1");
}

// ---------------------------------------------------------------------------
// Configuration file constructor
// ---------------------------------------------------------------------------

/// Type-erased pointer variant for generic parameter assignment.
using ParamPtr = std::variant<double*, long int*, int*>;

/**
 * @brief Reads parameter values from a key-value configuration file.
 *
 * Each line of the file should contain a key and a value separated by
 * whitespace.  Recognised keys are mapped to the corresponding member
 * variable.  Unknown, missing, duplicate, and invalid keys cause an exception.
 * Lines may contain # comments and an optional equals sign.
 *
 * After parsing, all derived quantities (box size, particle count, etc.)
 * are recomputed.
 *
 * @param path        Directory containing the configuration file.
 * @param configFile  Name of the configuration file.
 */
SimParams::SimParams(std::string path, std::string configFile) : SimParams() {
    // Map config-file keys to member-variable pointers
    std::map<std::string, ParamPtr> params = {
        {"rc",         &rc},
        {"sigma",      &sigma},
        {"epsilon",    &epsilon},
        {"mass1",      &mass1},
        {"mass2",      &mass2},
        {"NPartRatio", &NPartRatio},
        {"dt",         &dt},
        {"T",          &T},
        {"totsteps",   &totsteps},
        {"nx",         &nx},
        {"ny",         &ny},
        {"nz",         &nz},
        {"rho",        &rho},
    };

    asciiFileReader reader(path, configFile);
    std::vector<std::string> Lines = reader.getAllLines();

    std::set<std::string> seen;
    for (std::string line : Lines) {
        line = line.substr(0, line.find('#'));
        std::replace(line.begin(), line.end(), '=', ' ');
        std::istringstream iss(line);
        std::string key, value, extra;
        if (!(iss >> key)) continue;
        if (!(iss >> value) || (iss >> extra)) {
            throw std::invalid_argument("Expected key and value in config: " + line);
        }
        // Retained for compatibility with historical config files; CSV uses ','.
        if (key == "delimiter" && value == "0") continue;
        const auto it = params.find(key);
        if (it == params.end()) throw std::invalid_argument("Unknown config key: " + key);
        if (!seen.insert(key).second) throw std::invalid_argument("Duplicate config key: " + key);
        try {
            std::visit([&value](auto* arg) {
                using Value = std::decay_t<decltype(*arg)>;
                size_t used = 0;
                if constexpr (std::is_same_v<Value, double>) {
                    *arg = std::stod(value, &used);
                    if (!std::isfinite(*arg)) throw std::invalid_argument("non-finite value");
                } else if constexpr (std::is_same_v<Value, long int>) {
                    *arg = std::stol(value, &used);
                } else {
                    *arg = std::stoi(value, &used);
                }
                if (used != value.size()) throw std::invalid_argument("trailing characters");
            }, it->second);
        } catch (const std::exception&) {
            throw std::invalid_argument("Invalid config value for " + key + ": " + value);
        }
    }
    for (const auto& entry : params) {
        if (!seen.count(entry.first)) throw std::invalid_argument("Missing config key: " + entry.first);
    }
    if (nx <= 0 || ny <= 0 || nz <= 0 || rho <= 0 || sigma <= 0 ||
        epsilon <= 0 || rc <= 0 || mass1 <= 0 || mass2 <= 0 || dt <= 0 ||
        T < 0 || totsteps < 0 || NPartRatio < 0 || NPartRatio > 1) {
        throw std::invalid_argument("Invalid physical parameters: dimensions, density, LJ parameters, masses and dt must be positive; T and totsteps nonnegative; NPartRatio in [0,1].");
    }

    // Derive dependent quantities
    setNparticles();
    setNPartSmallMass();
    setTotTime();
    setRatio();
    setLx();
    setLy();
    setLz();

    if (!std::isfinite(Lx) || !std::isfinite(Ly) || !std::isfinite(Lz) ||
        !std::isfinite(totTime) || rc > 0.5 * std::min({Lx, Ly, Lz})) {
        throw std::invalid_argument("Require finite box/time and rc <= half the shortest box dimension.");
    }
    setOutputTimePath(path + "/OutputTimes");
    setEnergyPath(path);
}

SimParams::~SimParams() {
}

// ---------------------------------------------------------------------------
// Time-stepping parameter accessors
// ---------------------------------------------------------------------------

double SimParams::getDt() const { return dt; }
void SimParams::setDt(double dt) { this->dt = dt; }

long int SimParams::getTotsteps() const { return totsteps; }
void SimParams::setTotsteps(long int totsteps) { this->totsteps = totsteps; }

double SimParams::getTotTime() const { return totTime; }
void SimParams::setTotTime() { this->totTime = getTotsteps() * getDt(); }

// ---------------------------------------------------------------------------
// LJ interaction parameter accessors
// ---------------------------------------------------------------------------

double SimParams::getEpsilon() const { return epsilon; }
void SimParams::setEpsilon(double epsilon) { this->epsilon = epsilon; }

double SimParams::getSigma() const { return sigma; }
void SimParams::setSigma(double sigma) { this->sigma = sigma; }

double SimParams::getRc() const { return rc; }
void SimParams::setRc(double rc) { this->rc = rc; }

// ---------------------------------------------------------------------------
// Mass accessors
// ---------------------------------------------------------------------------

double SimParams::getMass1() const { return mass1; }
void SimParams::setMass1(double mass1) { this->mass1 = mass1; }

double SimParams::getMass2() const { return mass2; }
void SimParams::setMass2(double mass2) { this->mass2 = mass2; }

double SimParams::getNPartRatio() const { return NPartRatio; }
void SimParams::setNPartRatio(double nratio) { NPartRatio = nratio; }

// ---------------------------------------------------------------------------
// Thermodynamic parameter accessors
// ---------------------------------------------------------------------------

double SimParams::getT() const { return T; }
void SimParams::setT(double t) { this->T = t; }

double SimParams::getRho() const { return rho; }
void SimParams::setRho(double rho) { this->rho = rho; }

// ---------------------------------------------------------------------------
// Grid dimension accessors
// ---------------------------------------------------------------------------

int SimParams::getNx() const { return nx; }
void SimParams::setNx(int nx) { this->nx = nx; }

int SimParams::getNy() const { return ny; }
void SimParams::setNy(int ny) { this->ny = ny; }

int SimParams::getNz() const { return nz; }
void SimParams::setNz(int nz) { this->nz = nz; }

// ---------------------------------------------------------------------------
// Derived box dimensions
// ---------------------------------------------------------------------------

double SimParams::getLx() const { return Lx; }

/**
 * @brief Computes the x-dimension of the simulation box.
 *
 * The cell length is derived from the simple-cubic reference density
 * (rhomax = 1/sigma^3) and the target density rho:
 *   cell_length = sigma * (rhomax / rho)^(1/3)
 * Then Lx = nx * cell_length.
 */
void SimParams::setLx() {
    double vcell  = getSigma() * getSigma() * getSigma();
    double rhomax = 1.0 / vcell;
    double sf     = std::pow(rhomax / getRho(), 1.0 / 3);
    double clen   = sf * getSigma();
    Lx = getNx() * clen;
}

double SimParams::getLy() const { return Ly; }

void SimParams::setLy() {
    double vcell  = getSigma() * getSigma() * getSigma();
    double rhomax = 1.0 / vcell;
    double sf     = std::pow(rhomax / getRho(), 1.0 / 3);
    double clen   = sf * getSigma();
    Ly = getNy() * clen;
}

double SimParams::getLz() const { return Lz; }

void SimParams::setLz() {
    double vcell  = getSigma() * getSigma() * getSigma();
    double rhomax = 1.0 / vcell;
    double sf     = std::pow(rhomax / getRho(), 1.0 / 3);
    double clen   = sf * getSigma();
    Lz = getNz() * clen;
}

// ---------------------------------------------------------------------------
// Derived particle counts
// ---------------------------------------------------------------------------

long int SimParams::getNparticles() const { return Nparticles; }
void SimParams::setNparticles() {
    long int count = 1;
    for (int n : {getNx(), getNy(), getNz()}) {
        if (n <= 0 || count > std::numeric_limits<long int>::max() / n) {
            throw std::invalid_argument("Invalid or overflowing particle count.");
        }
        count *= n;
    }
    Nparticles = count;
}

long int SimParams::getNPartSmallMass() const { return NPartSmallMass; }
void SimParams::setNPartSmallMass() {
    NPartSmallMass = (long int)(getNparticles() * getNPartRatio());
}

int SimParams::getRatio() const { return Ratio; }

/**
 * @brief Computes the integer time-step ratio sqrt(m_heavy / m_light).
 *
 * Used by the multiple-time-step integrator to determine how many
 * sub-steps are taken for light particles per heavy-particle step.
 */
void SimParams::setRatio(void) {
    double m1 = getMass1();
    double m2 = getMass2();

    const double ratio = std::sqrt(std::max(m1, m2) / std::min(m1, m2));
    if (!std::isfinite(ratio) || ratio > std::numeric_limits<int>::max()) {
        throw std::invalid_argument("Mass ratio is too large.");
    }
    this->Ratio = std::max(1, static_cast<int>(ratio));
}

// ---------------------------------------------------------------------------
// Output path accessors
// ---------------------------------------------------------------------------

const std::string& SimParams::getOutputTimePath() const { return OutputTimePath; }
void SimParams::setOutputTimePath(const std::string& outputTimePath) {
    OutputTimePath = outputTimePath;
}

const std::string& SimParams::getEnergyPath() const { return EnergyPath; }
void SimParams::setEnergyPath(const std::string& energyPath) {
    EnergyPath = energyPath;
}

// ---------------------------------------------------------------------------
// Pretty-print
// ---------------------------------------------------------------------------

void SimParams::printSimParams(void) {
    printSimParamsPriv();
}

void SimParams::printSimParamsPriv(void) {
    std::cout << "Simulation Parameters:\n" << std::endl;
    std::cout << "nx             = " << nx << std::endl;
    std::cout << "ny             = " << ny << std::endl;
    std::cout << "nz             = " << nz << std::endl;
    std::cout << "Nparticles     = " << Nparticles << std::endl;
    std::cout << "NPartRatio     = " << NPartRatio << std::endl;
    std::cout << "NPartSmallMass = " << NPartSmallMass << std::endl;
    std::cout << "rc             = " << rc << std::endl;
    std::cout << "sigma          = " << sigma << std::endl;
    std::cout << "epsilon        = " << epsilon << std::endl;
    std::cout << "mass1          = " << mass1 << std::endl;
    std::cout << "mass2          = " << mass2 << std::endl;
    std::cout << "Ratio          = " << Ratio << std::endl;
    std::cout << "T              = " << T << std::endl;
    std::cout << "rho            = " << rho << std::endl;
    std::cout << "Lx             = " << Lx << std::endl;
    std::cout << "Ly             = " << Ly << std::endl;
    std::cout << "Lz             = " << Lz << std::endl;
    std::cout << "dt             = " << dt << std::endl;
    std::cout << "totsteps       = " << totsteps << std::endl;
    std::cout << "totTime        = " << totTime << std::endl;
    std::cout << "\n" << std::endl;
}
