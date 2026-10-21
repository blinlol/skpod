#include "solution.hpp"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void printUsage(const char* program) {
    std::cerr << "Usage: " << program << " N Lx Ly Lz T K\n";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 7) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    try {
        int n = std::stoi(argv[1]);
        double lx = std::stod(argv[2]);
        double ly = std::stod(argv[3]);
        double lz = std::stod(argv[4]);
        double t = std::stod(argv[5]);
        int k = std::stoi(argv[6]);

        Solution solution(lx, ly, lz, t, n, k);
        solution.run();

        std::cout << std::scientific << std::setprecision(8)
                  << "tau=" << (t / k)
                  << " time_sec=" << solution.seconds()
                  << " delta=" << solution.maxError() << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
