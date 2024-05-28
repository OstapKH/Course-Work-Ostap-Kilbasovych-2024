#ifndef AUTOMATAMPI_HPP_
#define AUTOMATAMPI_HPP_

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <limits>
#include <mpi.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>

#include "Cell.hpp"

#define D 2.0	// Euclidean spatial dimension (2 dimensions)

class AutomataMPI {
public:
// CGAL type definitions
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Delaunay_triangulation_2<K> Delaunay;
	typedef K::Point_2 CGAL_Point;
	
// Point type definition
	struct Point {
		double x, y, z;

		Point(double _x, double _y, double _z): x(_x), y(_y), z(_z){}
	};
	
// Tumor simulation parameters type definition
	struct TumorParameters {
		double tumorInitRadius;
		double a;				// Base necrotic thickness, controlled by nutritional needs (0.58 mm^(1/2))
		double b;				// Base proliferative thickness, controlled by nutritional needs (0.30 mm^(1/2))
		double p0;				// Base probability of division, linked to cell-doubling time (0.192 and 0.384)
		double gamma;			// Mutation rate (determines the number of invasive cells, 0.05)
		double xDegradation;	// ECM degradation ability (0.4 - 1.0)
		int muMotility;		// Cell motility (the number of "jumps" from one automaton cell to another, 1 - 3)
		
		TumorParameters(double tumorInitRadius=0.15, double a=0.58, double b=0.30, double p0=0.288, double gamma=0.05, double xDegradation=0.7, int muMotility=3):tumorInitRadius(tumorInitRadius),a(a),b(b),p0(p0),gamma(gamma),xDegradation(xDegradation),muMotility(muMotility){}
	};

// Proliferation cell and probability type definition
	struct ProliferationProbability {
		double pDiv;
		int ecmCellIndex;
	};

// Constructors
	AutomataMPI(int size, int rank, unsigned long int maxCells = 10000, int randSeed=1, double height=5.0, double width=5.0, double deep=5.0, double cellRadius=0.0075, TumorParameters tumorParam={0.15, 0.58});

// Destructor
	~AutomataMPI(){}
	
// Methods and functions
	double euclideanDistance(const Point& p1, const Point& p2);	// Calculates the distance between two points
	void storePointsToFile();
	void storeLastStateToFile();
	const std::vector<Cell>& getLastState(){return state1;}
	const std::vector<Cell>& computeNextState();
	
// Tumor metrics
	double asphericity();

private:
// Private methods
	void allocateCellsRandomly();
	void searchNeighbors();
	void updateTumorBorderCondition();
	
	bool quiescentToNectrotic(int cellIndex);		// Quiescent cells more than a certain distance deltaN from the tumor’s edge are turned necrotic.
	bool proliferativeToQuiescent(int cellIndex);	// Proliferatice cells more than a certain distance deltaP from the tumor’s edge are turned quiescent.
	ProliferationProbability proliferativeProliferates(int cellIndex, ProliferationProbability& res);	// Proliferates given proliferative cell to surrounding ECM. Returns the proliferate probability and the selected ECM cell index.
	int invasiveMoves(int cellIndex, int vectorIndex);	// Given invasive cell, it degrades the surrounding ECM and moves in search of nutrients and oxygen. vectorIndex is neccesary to update the value on tumorCellsIndex1.

// Attributes
	int size, rank;						// MPI attributes
	
	unsigned long int maxCells, maxAttempts;
	unsigned int nCells;				// Cell counters and parameters required for the space generation
	unsigned int nonInvasiveTumorCells;	// Non-invasive tumor cells counter
	unsigned int invasiveTumorCells;	// Invasive tumor cells counter
	Point centroid;						// Centroid of the tumoral cells
	int randSeed;						// Random seed for cell distribution
	double height, width, deep;			// Space size
	double cellRadius;					// Minimum radius of each cell
	TumorParameters tumorParam;			// Tumor simulation parameters
	std::vector<Cell> state0, state1;	// States of the automaton, 0 is the old state and 1 the next one
	Delaunay dt;						// Computes the neighborhood of each cell
	std::vector<int> tumorCellsIndex0;	// Registers the index of every tumoral cell for old state
	std::vector<int> tumorCellsIndex1;	// Registers the index of every tumoral cell for new state
	int stateIteration;					// Counts how many times the computeNextState method has been called
};

inline double AutomataMPI::euclideanDistance(const Point& p1, const Point& p2){
	return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

#endif /* AUTOMATAMPI_HPP_ */
