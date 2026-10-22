// Independent numerical and I/O regression checks, compiled against every variant.
#include "MDPart/MDPart.h"
#include "SimParams/SimParams.h"
#include "IO/ASCII/asciiFileWriter.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
void near(double actual, double expected, double tolerance, const char* message) {
    require(std::isfinite(actual) && std::abs(actual - expected) <= tolerance, message);
}
template<class F> void rejects(F operation, const char* message) {
    try { operation(); } catch (const std::exception&) { return; }
    throw std::runtime_error(message);
}
MDPart particle(long long id, double x) {
    MDPart p(id, 2.0, x, 0, 0, 0, 0, 0);
    p.setL(Eigen::Vector3d(10, 10, 10));
    p.setSigma(1.2); p.setEpsilon(0.8); p.setRc(2.5); p.setUcut();
    return p;
}
}

int main(int argc, char** argv) {
    try {
        require(argc == 2, "Need temporary output directory");
        const std::string path = argv[1];
        MDPart a = particle(0, 0), b = particle(1, 1.4);
        near(a.getAcceleration().norm(), 0, 0, "Default acceleration not zero");
        MDPart accelerated(9, 2.0, 0, 0, 0, 0, 0, 0, 1, 2, 3);
        near(accelerated.getAz(), 3, 0, "Constructor lost explicit acceleration");
        near(accelerated.getForceZ(), 6, 0, "Constructor force inconsistent with acceleration");
        const double q6 = std::pow(1.2 / 1.4, 6);
        const double qcut6 = std::pow(1.2 / 2.5, 6);
        const double expectedU = 4 * 0.8 * (q6*q6 - q6 - qcut6*qcut6 + qcut6);
        const double expectedFx = -24 * 0.8 * (2*q6*q6 - q6) / 1.4;
        double accumulated = 7.0;
        const auto force = a.computeForceijAndU(b, accumulated);
        near(force.x(), expectedFx, 1e-12, "Incorrect analytic LJ force");
        near(a.computeUij(b), expectedU, 1e-12, "Incorrect shifted potential");
        near(accumulated, 7.0 + expectedU, 1e-12, "Pair potential must accumulate");
        near((force + b.computeForceij(a)).norm(), 0, 1e-12, "Newton's third law failed");
        const double h = 1e-6;
        a.setX(h); const double plus = a.computeUij(b);
        a.setX(-h); const double minus = a.computeUij(b);
        a.setX(0);
        near(force.x(), -(plus-minus)/(2*h), 1e-7, "Force is not -gradient(U)");
        b.setX(11.4);
        near(a.computeUij(b), expectedU, 1e-12, "Minimum image failed");
        b.setX(2.5);
        near(a.computeUij(b), 0, 1e-12, "Potential not zero at cutoff");
        b.setX(2.6);
        near(a.computeForceij(b).norm(), 0, 0, "Force beyond cutoff");
        near(a.computeUij(b), 0, 0, "Potential beyond cutoff");
        a.setPosition(Eigen::Vector3d(21.3, -21.3, 10.2)); a.pbc();
        near(a.getX(), 1.3, 1e-12, "PBC x");
        near(a.getY(), -1.3, 1e-12, "PBC y");
        near(a.getZ(), 0.2, 1e-12, "PBC z");
        a.setForce(Eigen::Vector3d(2, 4, 6)); a.updateVel(0.5);
        near(a.getVz(), 1.5, 1e-12, "Kick acceleration F/m");
        rejects([] { MDPart p("0,1", ','); }, "Short CSV not rejected");
        rejects([] { MDPart p("0,1,0,0,0,0,0,nan", ','); }, "NaN CSV not rejected");
        rejects([] { MDPart p("0,1,0,0,0,0,0,1oops", ','); }, "CSV numeric suffix not rejected");
        rejects([] { PointInSpace p("0,1", ','); }, "Short point CSV not rejected");
        rejects([] { LagrangianPart p("0,1", ','); }, "Short Lagrangian CSV not rejected");
        SimParams defaults;
        require(std::isfinite(defaults.getLx()) && defaults.getNparticles() == 64,
                "Default parameters are invalid");
        {
            asciiFileWriter writer(path, "lines.txt");
            writer.writeOnFile(std::vector<std::string>{"zero", "one", "two", "three"});
        }
        asciiFileReader reader(path, "lines.txt");
        require(reader.getLineAt(1) == "one", "Seek row 1");
        require(reader.getLineAt(2) == "two", "Consecutive seek skips row");
        require(reader.getLineAt(2) == "two", "Repeated seek wrong");
        require(reader.getNextLine() == "three", "Sequential read after seek");
        require(reader.getLineAt(0) == "zero", "Backward seek");
        require(reader.getMultipleLines(1, 2) == std::vector<std::string>({"one", "two"}),
                "Range cursor wrong");
        require(reader.getAllLines().size() == 1, "Bulk read cursor wrong");
        require(reader.getNextLine() == "zero", "Bulk read must reset cursor");
        rejects([&] { reader.getLineAt(10); }, "Seek beyond EOF must fail");
        rejects([&] { asciiFileReader missing(path, "absent.txt"); }, "Missing file not rejected");
        {
            asciiFileWriter append(path, "lines.txt", std::ios::app);
            append.writeOnFile("four");
        }
        asciiFileReader afterAppend(path, "lines.txt");
        require(afterAppend.getAllLines().size() == 5, "Explicit append mode truncated file");
        std::cout << "Physics and I/O checks passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
