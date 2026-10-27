#include "solution.hpp"

#include <mpi.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {

void printUsage(const char* program) {
    std::cerr << "Usage: " << program << " N Lx Ly Lz T K\n";
}

}  // namespace

int main(int argc, char** argv) {
    int provided = 0;
    if (MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) != MPI_SUCCESS) {
        return EXIT_FAILURE;
    }

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int status = EXIT_SUCCESS;
    try {
        if (argc != 7) {
            if (rank == 0) {
                printUsage(argv[0]);
            }
            status = EXIT_FAILURE;
        } else {
            int n = std::stoi(argv[1]);
            double lx = std::stod(argv[2]);
            double ly = std::stod(argv[3]);
            double lz = std::stod(argv[4]);
            double t = std::stod(argv[5]);
            int k = std::stoi(argv[6]);

            Solution solution(lx, ly, lz, t, n, k);
            solution.run();
            double delta = solution.maxError();
            if (solution.rank() == 0) {
                std::cout << std::scientific << std::setprecision(8)
                          << "tau=" << (t / k)
                          << " time_sec=" << solution.seconds()
                          << " delta=" << delta << '\n';
            }
        }
    } catch (const std::exception& error) {
        if (rank == 0) {
            std::cerr << error.what() << '\n';
            std::cerr.flush();
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    MPI_Finalize();
    return status;
}
