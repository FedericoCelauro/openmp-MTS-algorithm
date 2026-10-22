/**
 * Trace the real MTS evolver and independently reconstruct energy each macro
 * step. State checkpoints expose the raw integrator buffers, without replacing
 * them with recomputed snapshot forces. No integration formulas are changed.
 */
#include "MDParticleEvolver/MDParticleEvolver.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>

int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::invalid_argument("Usage: trace <job-directory> <output-prefix>");
        SimParams params(argv[1], "config.txt");
        MDParticleEvolver evolver(params, argv[1], "initial_conf.csv", ',');
        auto& light = evolver.getParticlesMap().at(evolver.getLightMass());
        auto& heavy = evolver.getParticlesMap().at(evolver.getHeavyMass());
        std::vector<MDPart*> particles;
        for (auto& group : evolver.getParticlesMap()) {
            for (auto& particle : group.second) particles.push_back(&particle);
        }
        std::ofstream energy(std::string(argv[2])+".energy.csv");
        std::ofstream state(std::string(argv[2])+".state.csv");
        energy.exceptions(std::ios::failbit | std::ios::badbit);
        state.exceptions(std::ios::failbit | std::ios::badbit);
        energy << std::setprecision(17) << "step,K,U,K_fixed,U_fixed,E,E_fixed\n";
        state << std::setprecision(17)
              << "step,id,m,x,y,z,vx,vy,vz,fx,fy,fz,fxAux,fyAux,fzAux\n";
        const std::set<long> checkpoints = {0,1,2,5,10,25,50,100,200,300,400,600,800,1000};
        for (long step = 0; step <= params.getTotsteps(); ++step) {
            long double kinetic = 0.0L, potential = 0.0L;
            for (size_t i = 0; i < particles.size(); ++i) {
                const auto& p = *particles[i];
                kinetic += 0.5L*p.getMass()*p.getVelocity().squaredNorm();
                for (size_t j = i+1; j < particles.size(); ++j) {
                    double pairEnergy = 0.0;
                    particles[i]->computeForceijAndU(*particles[j], pairEnergy);
                    potential += pairEnergy;
                }
            }
            const double reportedK = evolver.computeKineticEnergy();
            const double reportedU = evolver.getU();
            energy << step << ',' << reportedK << ',' << reportedU << ','
                   << kinetic << ',' << potential << ',' << reportedK+reportedU << ','
                   << kinetic+potential << '\n';
            if (checkpoints.count(step) || step == params.getTotsteps()) {
                for (const auto* p : particles) {
                    state << step << ',' << p->getId() << ',' << p->getMass() << ','
                          << p->getX() << ',' << p->getY() << ',' << p->getZ() << ','
                          << p->getVx() << ',' << p->getVy() << ',' << p->getVz() << ','
                          << p->getForceX() << ',' << p->getForceY() << ',' << p->getForceZ() << ','
                          << p->getForceXAux() << ',' << p->getForceYAux() << ',' << p->getForceZAux() << '\n';
                }
            }
            if (step == params.getTotsteps()) break;
            evolver.evolve(light, heavy, params.getDt());
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
