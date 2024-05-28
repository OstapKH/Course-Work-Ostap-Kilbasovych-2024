#include <iostream>
#include <cstdlib>

#include "Automata.hpp"

using namespace std;

int main(int argc, char **argv){
	int rank, size;

	if(rank==0){
		unsigned long int maxIter = 1;
		unsigned long int maxCells = 10000;
		int seed = 42;
		double size = 5.0;			// milimeters
		double cellRadius = 0.0075;	// milimeters
		double tumorRadius = 0.15;	// milimeters
		double a = 0.58;			// milimeters^(1/2)
		double b = 0.30;			// milimeters^(1/2)
		double p0 = 0.288;			// probability
		double gamma = 0.05;		// probability
		double xDegradation = 0.7;	// probability
		int muMotility = 3;		// jumps

		if(argc > 1){
			maxIter = strtoul(argv[1], nullptr, 10);
			if(argc > 2){
				maxCells = strtoul(argv[2], nullptr, 10);
				if(argc > 3){
					seed = atoi(argv[3]);
					if(argc > 4){
						size = strtod(argv[4], nullptr);
						if(argc > 5){
							cellRadius = strtod(argv[5], nullptr);
							if(argc > 6){
								tumorRadius = strtod(argv[6], nullptr);
								if(argc > 7){
									a = strtod(argv[7], nullptr);
									if(argc > 8){
										b = strtod(argv[8], nullptr);
										if(argc > 9){
											p0 = strtod(argv[9], nullptr);
											if(argc > 10){
												gamma = strtod(argv[10], nullptr);
												if(argc > 11){
													xDegradation = strtod(argv[11], nullptr);
													if(argc > 12){
														muMotility = atoi(argv[12]);
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}

		cout
		<<"Arguments:\n\t- Maximum number of iterations: "<<maxIter
		<<"\n\t- Maximum number of cells: "<<maxCells
		<<"\n\t- Random seed: "<<seed
		<<"\n\t- Space size (mm): "<<size
		<<"\n\t- Minimum cell radius (mm): "<<cellRadius
		<<"\n\t- Initial tumor radius (mm): "<<tumorRadius
		<<"\n\t- Constant \"a\" - necrotic thickness (mm^(1/2)): "<<a
		<<"\n\t- Constant \"b\" - proliferative thickness (mm^(1/2)): "<<b
		<<"\n\t- Constant \"p0\" - division rate (probability): "<<p0
		<<"\n\t- Constant \"gamma\" - mutation rate (probability): "<<gamma
		<<"\n\t- Constant \"X\" - ECM degradation ability (probability): "<<xDegradation
		<<"\n\t- Constant \"mu\" - cell motility (jumps): "<<muMotility
		<<endl;

		Automata automata(maxCells, seed, size, size, size, cellRadius, {tumorRadius,a,b,p0,gamma,xDegradation,muMotility});
		automata.storeLastStateToFile();

		for(unsigned long int i=0; i<maxIter; i++){
			automata.computeNextState();
			if(i%10 == 0){
				automata.storeLastStateToFile();
			}
		}
	}
	return 0;
}

