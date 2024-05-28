#define INCLUDE_MPI

#include <iostream>
#include <cstdlib>
#include "DiffusionReaction3D.hpp"

#ifdef INCLUDE_MPI

#include <mpi.h>

#endif

using namespace std;

int main(int argc, char **argv) {
#ifdef INCLUDE_MPI
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
#endif
        int domainSize = 10;
        int numIterations = 100;
        bool isVascularized = false;

        for (int i = 1; i < argc; i++) {
            switch (i) {
                case 1:
                    domainSize = strtoul(argv[i], nullptr, 10);
                    break;
                case 2:
                    numIterations = strtoul(argv[i], nullptr, 10);
                    break;
                case 3:
                    isVascularized = static_cast<bool>(atoi(argv[i]));
                    break;
                default:
                    break;
            }
        }

        cout << "Arguments:\n\t- Domain size: " << domainSize
             << "\n\t- Number of iterations: " << numIterations
             << "\n\t- Is vascularized: " << isVascularized
             << endl;

#ifdef INCLUDE_MPI
        DiffusionReaction3D diffusionReaction3D(rank, size, domainSize, numIterations, isVascularized);
#else
        DiffusionReaction3D diffusionReaction3D(domainSize, numIterations, isVascularized);
#endif
        diffusionReaction3D.storeLastStateToFile();

        for (unsigned long int i = 0; i < numIterations; i++) {
            diffusionReaction3D.computeNextState();
        }


#ifdef INCLUDE_MPI
    } else {
        DiffusionReaction3D automata(size, true); // Placeholder for non-root processes
    }
    MPI_Finalize();
#endif

    return 0;
}
