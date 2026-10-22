/**
 * Time only calls to the real particle evolver. Configuration parsing,
 * initial-force calculation, energy checks and output are outside the timer.
 * Build against each variant's existing releaseO3 objects, excluding main.
 */
#include "MDParticleEvolver/MDParticleEvolver.h"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::invalid_argument("Usage: benchmark <job-directory>");
        SimParams params(argv[1], "config.txt");
        MDParticleEvolver evolver(params, argv[1], "initial_conf.csv", ',');
        const double initial = evolver.computeKineticEnergy() + evolver.getU();
        int threads = 1;
#ifdef _OPENMP
        // Start workers before timing and record the actual team size.
        #pragma omp parallel
        {
            #pragma omp single
            threads = omp_get_num_threads();
        }
#endif
#ifdef MTS
        auto& light = evolver.getParticlesMap().at(evolver.getLightMass());
        auto& heavy = evolver.getParticlesMap().at(evolver.getHeavyMass());
#endif
        const auto start = std::chrono::steady_clock::now();
        for (long int step = 0; step < params.getTotsteps(); ++step) {
#ifdef MTS
            evolver.evolve(light, heavy, params.getDt());
#else
            evolver.evolve(params.getDt());
#endif
        }
        const auto stop = std::chrono::steady_clock::now();
        const double final = evolver.computeKineticEnergy() + evolver.getU();
        if (!std::isfinite(initial) || !std::isfinite(final)) {
            throw std::runtime_error("Non-finite energy in benchmark");
        }
        std::cout << std::setprecision(17)
                  << "{\"seconds\":" << std::chrono::duration<double>(stop-start).count()
                  << ",\"threads\":" << threads
                  << ",\"initial_energy\":" << initial << ",\"final_energy\":" << final << "}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
