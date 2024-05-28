#define INCLUDE_MPI

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <array>
#include <unordered_set>
#include <cstdlib>

#ifdef INCLUDE_MPI
#include <mpi.h>
#include "AutomataMPI3D.hpp"
#else
#include "Automata.hpp"
#endif



using namespace std;

int main(int argc, char **argv){
#ifdef INCLUDE_MPI
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if(rank==0){
#endif
		unsigned long int maxIter = 1;
		unsigned long int maxCells = 10000;
		int seed = 42;
		double spaceSize = 5.0;		// milimeters
		double cellRadius = 0.0075;	// milimeters
		double tumorRadius = 0.15;	// milimeters
		double a = 0.58;			// milimeters^(1/2)
		double b = 0.30;			// milimeters^(1/2)
		double p0 = 0.288;			// probability
		double gamma = 0.05;		// probability
		double xDegradation = 0.7;	// probability
		int muMotility = 3;			// jumps
		
		for(int i=1; i< argc; i++){
			switch(i){
			case 1:
				maxIter = strtoul(argv[i], nullptr, 10);
				break;
			case 2:
				maxCells = strtoul(argv[i], nullptr, 10);
				break;
			case 3:
				seed = atoi(argv[i]);
				break;
			case 4:
				spaceSize = strtod(argv[i], nullptr);
				break;
			case 5:
				cellRadius = strtod(argv[i], nullptr);
				break;
			case 6:
				tumorRadius = strtod(argv[i], nullptr);
				break;
			case 7:
				a = strtod(argv[i], nullptr);
				break;
			case 8:
				b = strtod(argv[i], nullptr);
				break;
			case 9:
				p0 = strtod(argv[i], nullptr);
				break;
			case 10:
				gamma = strtod(argv[i], nullptr);
				break;
			case 11:
				xDegradation = strtod(argv[i], nullptr);
				break;
			case 12:
				muMotility = atoi(argv[i]);
				break;
				
			default:
				break;
			}
		}
		
		cout
		<<"Arguments:\n\t- Maximum number of iterations: "<<maxIter
		<<"\n\t- Maximum number of cells: "<<maxCells
		<<"\n\t- Random seed: "<<seed
		<<"\n\t- Space size (mm): "<<spaceSize
		<<"\n\t- Minimum cell radius (mm): "<<cellRadius
		<<"\n\t- Initial tumor radius (mm): "<<tumorRadius
		<<"\n\t- Constant \"a\" - necrotic thickness (mm^(1/2)): "<<a
		<<"\n\t- Constant \"b\" - proliferative thickness (mm^(1/2)): "<<b
		<<"\n\t- Constant \"p0\" - division rate (probability): "<<p0
		<<"\n\t- Constant \"gamma\" - mutation rate (probability): "<<gamma
		<<"\n\t- Constant \"X\" - ECM degradation ability (probability): "<<xDegradation
		<<"\n\t- Constant \"mu\" - cell motility (jumps): "<<muMotility
		<<endl;

#ifdef INCLUDE_MPI		
		AutomataMPI3D automata(size, rank, maxCells, seed, spaceSize, spaceSize, spaceSize, cellRadius, {tumorRadius,a,b,p0,gamma,xDegradation,muMotility});
#else
		Automata automata(maxCells, seed, spaceSize, spaceSize, spaceSize, cellRadius, {tumorRadius,a,b,p0,gamma,xDegradation,muMotility});
#endif
		automata.storeLastStateToFile();
		
		for(unsigned long int i=0; i<maxIter; i++){
			automata.computeNextState();
			if(i%10 == 0){
				automata.storeLastStateToFile();
			}
		}

#ifdef INCLUDE_MPI
	}else{
		AutomataMPI3D automata(size, rank);
	}
	MPI_Finalize();
#endif

	return 0;
}

