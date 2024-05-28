#include <iostream>
#include <fstream>
#include <cstdlib>
#include <mpi.h>

#define GRID_SIZE 100
#define BOUNDARY_VALUE 1.0

using namespace std;

int main(int argc, char **argv){
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	const int width=GRID_SIZE, height=GRID_SIZE;
	double fi0, D0, lambda14, Dt, Dx, Kmet, nCA;
	double **domain = new double*[height];
	double **domainNext = new double*[height];
	double **auxPtr;
	for(int i=0; i<height; i++){
		domain[i] = new double[width];
		domainNext[i] = new double[width];
	}

	fi0 = 1;
	D0 = 1.3e-6;
	lambda14 = 2.5e-4;
	Dt = 3600;
	Dx = 0.01;
	Kmet = 2.0e-4;
	nCA = 1;

	if(rank==0){
		for(int i=0; i<height; i++){
			for(int j=0; j<width; j++){
				domain[i][j] = (i==0 || i==height-1 || j==0 || j==width-1)? BOUNDARY_VALUE : 0;	// If it's a boundary cell, is initialized as 1, which means that has drugs, in other case 0
			}
		}

		for(int day=0; day<120; day++){
			ofstream file("day"+(day<10?string("0"):string(""))+to_string(day)+".txt");
			for(int i=0; i<height; i++){
				for(int j=0; j<width; j++){
					domainNext[i][j] = domain[i][j] + Dt*((D0*((i+1>=height ? BOUNDARY_VALUE : domain[i+1][j]) + (i-1<0 ? BOUNDARY_VALUE : domain[i-1][j]) + (j+1>=width ? BOUNDARY_VALUE : domain[i][j+1]) + (j-1<0 ? BOUNDARY_VALUE : domain[i][j-1]) - (4*domain[i][j]))/(Dx*Dx)) - Kmet*domain[i][j] - lambda14*nCA*(domain[i][j]/(domain[i][j] + fi0)));
					file<<domainNext[i][j]<<" ";
				}
				file<<endl;
			}
			file.close();
			auxPtr = domain;
			domain = domainNext;
			domainNext = auxPtr;
		}
	}

	MPI_Finalize();
	return 0;
}
