/**
 * @file MolecularDynamics.cpp
 * @brief Implementation of the MolecularDynamics class and program entry point.
 *
 * This file contains:
 *   - The main() function that parses command-line arguments, creates the
 *     MolecularDynamics object, runs the simulation, and reports wall-clock time.
 *   - The MolecularDynamics class methods for initialization and execution.
 *
 * @author Federico Celauro
 * @date   Apr 8, 2026
 */

#include "include/MolecularDynamics.h"

// ---------------------------------------------------------------------------
// Constructors & Destructor
// ---------------------------------------------------------------------------

MolecularDynamics::MolecularDynamics() {
    SimName = "sim1";
    pathInitialConditions = "jobs/sim1";

    if (!init()) {
        std::cerr << "Molecular Dynamics initialization failed!" << std::endl;
    }
}

MolecularDynamics::MolecularDynamics(char** argv) {
    SimName = argv[1];
    const fs::path requested(SimName);
    // Explicit paths are relative to the caller; job names resolve beside bin/.
    pathInitialConditions = fs::is_directory(requested)
        ? requested.string()
        : (fs::absolute(fs::path(argv[0])).parent_path().parent_path()
           / "jobs" / requested).string();

    if (!init()) {
        std::cerr << "Molecular Dynamics initialization failed!" << std::endl;
    }
}

MolecularDynamics::~MolecularDynamics() {
}

// ---------------------------------------------------------------------------
// Program entry point
// ---------------------------------------------------------------------------

/**
 * @brief Application entry point.
 *
 * Usage: ./MolecularDynamics <SimulationName>
 *
 * The simulation name corresponds to a subdirectory under jobs/ (beside bin/) that
 * must contain config.txt and initial_conf.csv.
 */
int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <SimulationName>" << std::endl;
        return 1;
    }

    auto start = std::chrono::high_resolution_clock::now();
    int exitCode = 0;

    try {
        MolecularDynamics simulation(argv);
        if (simulation.start()) {
            std::cout << "All done" << std::endl;
        } else {
            std::cerr << "Execution failed!" << std::endl;
            exitCode = 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Simulation failed: " << e.what() << std::endl;
        exitCode = 1;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "-------------------------------" << std::endl;
    std::cout << "Execution time: " << elapsed.count() << " seconds"
              << std::endl;
    std::cout << "-------------------------------" << std::endl;

    return exitCode;
}

// ---------------------------------------------------------------------------
// Component accessors
// ---------------------------------------------------------------------------

MDParticleEvolver* MolecularDynamics::getParticleEvolver() const {
    return particleEvolver.get();
}

void MolecularDynamics::setParticleEvolver(std::unique_ptr<MDParticleEvolver> evolver) {
    this->particleEvolver = std::move(evolver);
}

MDTimeIntegrator* MolecularDynamics::getTimeIntegrator() const {
    return timeIntegrator.get();
}

void MolecularDynamics::setTimeIntegrator(std::unique_ptr<MDTimeIntegrator> integrator) {
    this->timeIntegrator = std::move(integrator);
}

// ---------------------------------------------------------------------------
// Path accessors
// ---------------------------------------------------------------------------

const std::string& MolecularDynamics::getPathInitialConditions() const {
    return pathInitialConditions;
}

void MolecularDynamics::setPathInitialConditions(const std::string& pathInitialConditions) {
    this->pathInitialConditions = pathInitialConditions;
}

const std::string& MolecularDynamics::getConfigsFileName() const {
    return configsFileName;
}

void MolecularDynamics::setConfigsFileName(const std::string& configsFileName) {
    this->configsFileName = configsFileName;
}

const std::string& MolecularDynamics::getInitialConditionsFileName() const {
    return initialConditionsFileName;
}

void MolecularDynamics::setInitialConditionsFileName(const std::string& initialConditionsFileName) {
    this->initialConditionsFileName = initialConditionsFileName;
}

// ---------------------------------------------------------------------------
// Simulation name accessors
// ---------------------------------------------------------------------------

const std::string& MolecularDynamics::getSimName() const {
    return SimName;
}

void MolecularDynamics::setSimName(const std::string& simName) {
    SimName = simName;
}

// ---------------------------------------------------------------------------
// Simulation execution
// ---------------------------------------------------------------------------

bool MolecularDynamics::start() {
    bool result = timeIntegrator->run();
    return result;
}

// ---------------------------------------------------------------------------
// Initialization
// ---------------------------------------------------------------------------

/**
 * @brief Initializes all simulation components.
 *
 * Steps:
 *   1. Read simulation parameters from the config file.
 *   2. Print the parsed parameters for verification.
 *   3. Create the particle evolver (reads initial conditions).
 *   4. Create the time integrator.
 *   5. Write the energy CSV header.
 *
 * @return true if all components were initialized successfully.
 */
bool MolecularDynamics::init() {
    bool result = true;

    // Parse simulation parameters from config file
    SimParams params = SimParams(pathInitialConditions, configsFileName);
    params.printSimParams();

    // Initialize particle evolver from CSV initial conditions
    particleEvolver = std::make_unique<MDParticleEvolver>(
        params, pathInitialConditions, initialConditionsFileName, ',');

    // Initialize time integrator
    timeIntegrator = std::make_unique<MDTimeIntegrator>(
        params, particleEvolver.get());

    if (!particleEvolver) {
        std::cerr << "Failed to initialize MDParticleEvolver" << std::endl;
        result = false;
    }
    if (!timeIntegrator) {
        std::cerr << "Failed to initialize MDTimeIntegrator" << std::endl;
        result = false;
    }

    fs::create_directories(params.getOutputTimePath());

    // Write the energy CSV header
    asciiFileWriter fileWriter(params.getEnergyPath(), "energy.csv");
    char delimiter = ',';
    fileWriter.getPointer() << std::string("t") + delimiter
                            + "K" + delimiter
                            + "U" + delimiter
                            + "Etot" << std::endl;

    return result;
}
