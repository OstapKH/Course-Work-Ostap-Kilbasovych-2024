#include "Automata.hpp"

// ####################################################################################################################################################################
// Constructor |#######################################################################################################################################################
// ####################################################################################################################################################################
Automata::Automata(unsigned long int maxCells, int randSeed, double height, double width, double deep, double cellRadius, TumorParameters tumorParam, const std::vector<Cell> &cells):maxCells(maxCells),nCells(0),maxAttempts(10*maxCells),centroid(Point(0,0,0)),randSeed(randSeed),height(height),width(width),deep(deep),cellRadius(cellRadius),tumorParam(tumorParam),state0(std::vector<Cell>(cells)),state1(std::vector<Cell>(cells)){
	stateIteration = 0;
	dt = Delaunay();

	srand(randSeed);
	allocateCellsRandomly();
	searchNeighbors();
	updateTumorBorderCondition();
	
	nonInvasiveTumorCells = tumorCellsIndex1.size();
	invasiveTumorCells = 0;
	
	std::cout<<"Initial Tumoral Cells: "<<tumorCellsIndex1.size()<<std::endl;
	std::cout<<"Total number of Cells: "<<state1.size()<<std::endl;
}

// ####################################################################################################################################################################
// Finite States Machine of the automaton |############################################################################################################################
// ####################################################################################################################################################################
const std::vector<Cell>& Automata::computeNextState(){
	stateIteration++;
	
	state0 = state1;
	tumorCellsIndex0 = tumorCellsIndex1;
	
	//for(unsigned int i=0; i<nCells; i++){
	for(unsigned int j=0; j<tumorCellsIndex0.size(); j++){
		unsigned int i = tumorCellsIndex0[j];

		switch(state0[i].getType()){
		case Cell::NECROTIC:
		case Cell::ECM:
			break;
		case Cell::PROLIFERATIVE:
			ProliferationProbability res;
			proliferativeProliferates(i,res);
			if((static_cast<double>(rand())/RAND_MAX) < res.pDiv){
				state1[res.ecmCellIndex].setPECM(0.0);
				tumorCellsIndex1.push_back(res.ecmCellIndex);	// Registers the index of every new tumoral cell
				if((static_cast<double>(rand())/RAND_MAX) > tumorParam.gamma){
					state1[res.ecmCellIndex].setType(Cell::PROLIFERATIVE);	// Update of the condition of the new proliferative cell
					nonInvasiveTumorCells++;
					
					// Centroid Update
					centroid.x = (centroid.x * (nonInvasiveTumorCells-1)) + state1[res.ecmCellIndex].getX();
					centroid.y = (centroid.y * (nonInvasiveTumorCells-1)) + state1[res.ecmCellIndex].getY();
					//centroid.z = (centroid.z * (nonInvasiveTumorCells-1)) + state1[res.ecmCellIndex].getZ();
					centroid.x /= nonInvasiveTumorCells;
					centroid.y /= nonInvasiveTumorCells;
					//centroid.z /= nonInvasiveTumorCells;
				}else{
					state1[res.ecmCellIndex].setType(Cell::INVASIVE);	// Update of the condition of the new invasive cell
					invasiveTumorCells++;
				}
			}else if(Automata::proliferativeToQuiescent(i)){
				state1[i].setType(Cell::QUIESCENT);
			}
			break;
		case Cell::QUIESCENT:
			if(quiescentToNectrotic(i))
				state1[i].setType(Cell::NECROTIC);
			break;
		case Cell::INVASIVE:
			invasiveMoves(i,j);
			break;
		}
	}
	updateTumorBorderCondition();
	return state1;
}

// ####################################################################################################################################################################
// Next state functions and methods |##################################################################################################################################
// ####################################################################################################################################################################
bool Automata::quiescentToNectrotic(int cellIndex){
	int nearestBoundaryIndex = 0;								// Index of the nearest boundary cell
	double minDistance = std::numeric_limits<double>::max();	// Minimum distance from cell cellIndex and the boundary
	double Lt, deltaN;
	
	for(int n : tumorCellsIndex0){
		if(state0[n].getIsBoundary() && (euclideanDistance({state0[cellIndex].getX(),state0[cellIndex].getY(),0},{state0[n].getX(),state0[n].getY(),0}) < minDistance)){
			minDistance = euclideanDistance({state0[cellIndex].getX(),state0[cellIndex].getY(),0},{state0[n].getX(),state0[n].getY(),0});
			nearestBoundaryIndex = n;
		}
	}
	
	Lt = euclideanDistance({state0[nearestBoundaryIndex].getX(),state0[nearestBoundaryIndex].getY(),0}, centroid);
	
	deltaN = tumorParam.a * (pow(Lt, (double)(D-1)/D));
	
	return (minDistance > deltaN);
}

bool Automata::proliferativeToQuiescent(int cellIndex){
	int nearestBoundaryIndex = 0;								// Index of the nearest boundary cell
	double minDistance = std::numeric_limits<double>::max();	// Minimum distance from cell cellIndex and the boundary
	double Lt, deltaP;
	
	for(int n : tumorCellsIndex0){
		if(state0[n].getIsBoundary() && (euclideanDistance({state0[cellIndex].getX(),state0[cellIndex].getY(),0},{state0[n].getX(),state0[n].getY(),0}) < minDistance)){
			minDistance = euclideanDistance({state0[cellIndex].getX(),state0[cellIndex].getY(),0},{state0[n].getX(),state0[n].getY(),0});
			nearestBoundaryIndex = n;
		}
	}
	
	Lt = euclideanDistance({state0[nearestBoundaryIndex].getX(),state0[nearestBoundaryIndex].getY(),0}, centroid);
	
	deltaP = tumorParam.b * (pow(Lt, (double)(D-1)/D));
	
	return (minDistance > deltaP);
}

Automata::ProliferationProbability Automata::proliferativeProliferates(int cellIndex, ProliferationProbability& res){
	//ProliferationProbability res;
	res.pDiv = 0;										// Division probability
	res.ecmCellIndex = -1;								// Index of the nearest boundary cell
	double Lmax = std::numeric_limits<double>::max();	// Distance between the closest growth-permitting boundary cell in the direction of tumor growth and the tumor's geometric centroid
	double r = euclideanDistance({state0[cellIndex].getX(),state0[cellIndex].getY(),0}, centroid);	// Distance between the proliferative cell and the tumor centroid
	
	if(state0[cellIndex].getIsBoundary()){				// If isn't a tumor boundary cell, cannot proliferate, because doesn't have any ECM neighbor cell
		for(int n : state0[cellIndex].getNeighbors()){
			if((state0[n].getType() == Cell::ECM) && (euclideanDistance({state0[n].getX(),state0[n].getY(),0},centroid) < Lmax)){
				Lmax = euclideanDistance({state0[n].getX(),state0[n].getY(),0},centroid);
				res.ecmCellIndex = n;
			}
		}
		
		res.pDiv = (tumorParam.p0/2)*((1-(r/Lmax))+(1-state0[res.ecmCellIndex].getPECM()));
	}
	
	return res;
}

int Automata::invasiveMoves(int cellIndex, int vectorIndex){
	double deltaRho;	
	double bestECMCellDistance;
	int bestECMCellIndex = cellIndex;
	
	for(int m = 1 + (static_cast<int>(rand()) % tumorParam.muMotility); m > 0; m--){
		deltaRho = tumorParam.xDegradation * (static_cast<double>(rand())/RAND_MAX);
		for(int n : state0[cellIndex].getNeighbors()){
			if(state0[n].getType() == Cell::ECM){
				state1[n].setPECM(state0[n].getPECM() - deltaRho);
			}
		}
		
		bestECMCellDistance = std::numeric_limits<double>::min();
		for(int n : state1[cellIndex].getNeighbors()){
			if((state0[n].getType() == Cell::ECM) && (state1[n].getPECM() == 0.0) && (euclideanDistance({state0[n].getX(),state0[n].getY(),0},centroid) > bestECMCellDistance)){
				bestECMCellDistance = euclideanDistance({state0[n].getX(),state0[n].getY(),0},centroid);
				bestECMCellIndex = n;
			}
		}
		if(bestECMCellIndex != cellIndex){
			state1[cellIndex].setType(Cell::ECM);
			state1[bestECMCellIndex].setType(Cell::INVASIVE);
			cellIndex = bestECMCellIndex;
			tumorCellsIndex1[vectorIndex] = cellIndex;
		}
	}
	
	return 0;
}

// ####################################################################################################################################################################
// Tumor metrics methods |#############################################################################################################################################
// ####################################################################################################################################################################

double Automata::asphericity(){
	double inner = std::numeric_limits<double>::max();
	double outer = std::numeric_limits<double>::min();
	int innerCell, outerCell;
	
	for(int n: tumorCellsIndex1){
		if(state1[n].getIsBoundary() && (euclideanDistance({state1[n].getX(),state1[n].getY(),0},centroid) < inner)){
			inner = euclideanDistance({state1[n].getX(),state1[n].getY(),0},centroid);
			innerCell = n;
		}
		if(state1[n].getIsBoundary() && (euclideanDistance({state1[n].getX(),state1[n].getY(),0},centroid) > outer)){
			outer = euclideanDistance({state1[n].getX(),state1[n].getY(),0},centroid);
			outerCell = n;
		}
	}
	
	return outer/inner;
}

// ####################################################################################################################################################################
// Auxiliar and creation/initialization methods |######################################################################################################################
// ####################################################################################################################################################################
void Automata::allocateCellsRandomly(){
	// Generates random points until filling the space
	unsigned long int attempts = 0;
	Point center(width/2, height/2, deep/2);
	double x, y;
	bool matchesCondition;
	while(attempts < maxAttempts){
		// Generates random point inside the space
		x = static_cast<double>(rand()) / RAND_MAX * width;
		y = static_cast<double>(rand()) / RAND_MAX * height;

		// Verifies if the point is at a distance greater or equal to it's
		// radius with respect to all existing points
		matchesCondition = all_of(state1.begin(), state1.end(),
			[&](const Cell& cell) {
				return euclideanDistance({x,y,0}, {cell.getX(),cell.getY(),0}) >= 2 * cellRadius;
			});

		// Adds the point if matches the condition
		if(matchesCondition){
			if(euclideanDistance(center, {x,y,0}) <= tumorParam.tumorInitRadius){
				state1.emplace_back(x, y, 0, 0.0, false, Cell::PROLIFERATIVE);	// TODO: Initialize cells properly
				state0.emplace_back(x, y, 0, 0.0, false, Cell::PROLIFERATIVE);	// Both states must be initialized with the same data
				centroid.x += x;
				centroid.y += y;
				tumorCellsIndex1.push_back(nCells);	// Registers the index of every tumoral cell
				tumorCellsIndex0.push_back(nCells);	// Registers the index of every tumoral cell
			}else{
				state1.emplace_back(x, y, 0, 0.45);	// TODO: Initialize cells properly
				state0.emplace_back(x, y, 0, 0.45);	// Both states must be initialized with the same data
			}
			nCells++;
		}

		// Exits the loop if the desired cells number is reached
		if(state1.size() >= maxCells){
			break;
		}

		// Increase attempts counter
		attempts++;
	}
	
	centroid.x /= tumorCellsIndex1.size();
	centroid.y /= tumorCellsIndex1.size();
}

void Automata::searchNeighbors(){
	// Perform Delaunay triangulation with CGAL
	std::vector<CGAL_Point> points;
	for(const Cell& cell : state1){
		//points.push_back(CGAL_Point(cell.getX(), cell.getY()));
		points.emplace_back(CGAL_Point(cell.getX(), cell.getY()));
	}

	dt.insert(points.begin(), points.end());
	
	// Find neighbors list for each point
	int cellIndex, neighbor;
	for(Delaunay::Finite_vertices_iterator v = dt.finite_vertices_begin(); v != dt.finite_vertices_end(); ++v){
		cellIndex = std::distance(points.begin(), std::find(points.begin(), points.end(), v->point()));	// That's neccesary to find the position in the vector from the iterator
		Delaunay::Vertex_circulator c = dt.incident_vertices(v);
		if(c != 0){
			do{
				neighbor = std::distance(points.begin(), std::find(points.begin(), points.end(), c->point()));
					state1[cellIndex].getNeighbors().push_back(neighbor);	// That's neccesary to find the position in the vector from the iterator
			}while(++c != dt.incident_vertices(v));
		}
	}
}

void Automata::updateTumorBorderCondition(){
	for(int n : tumorCellsIndex1){
		if(state1[n].getType() != Cell::INVASIVE){
			for(int i : state1[n].getNeighbors()){	// Counts only the tumoral neighbors of the tumoral cell n
				if(state1[i].getType() == Cell::ECM){
					state1[n].setIsBoundary(true);
					break;
				}else{
					state1[n].setIsBoundary(false);
				}
			}
		}
	}
}

// ####################################################################################################################################################################
// Write files methods |###############################################################################################################################################
// ####################################################################################################################################################################
void Automata::storePointsToFile(){
	std::ofstream file("space.dat");
	for (Delaunay::Finite_vertices_iterator v = dt.finite_vertices_begin(); v != dt.finite_vertices_end(); ++v) {
		file << v->point().x() << " " << v->point().y() << std::endl;
	}
	file << "e" << std::endl;

	for (Delaunay::Finite_faces_iterator f = dt.finite_faces_begin(); f != dt.finite_faces_end(); ++f) {
		file << f->vertex(0)->point().x() << " " << f->vertex(0)->point().y() << std::endl;
		file << f->vertex(1)->point().x() << " " << f->vertex(1)->point().y() << std::endl;
		file << f->vertex(2)->point().x() << " " << f->vertex(2)->point().y() << std::endl;
		file << f->vertex(0)->point().x() << " " << f->vertex(0)->point().y() << std::endl;
		file << std::endl;
	}
	file.close();
	
#ifdef ENABLE_DRAWS
	CGAL::draw(dt);
#endif

	std::cout<<"Command to plot with Gnuplot: plot 'space.dat' with points title 'Points', '' with lines title 'Delaunay Triangulation'"<<std::endl;
}

void Automata::storeLastStateToFile(){
	std::ofstream file("state"+std::to_string(stateIteration)+".dat");
	for(Cell c : state1){
		file<<c.getX()<<" "<<c.getY()<<" ";
		switch(c.getType()){
		case Cell::ECM:
			file<<c.getPECM()<<std::endl;
			break;
		case Cell::PROLIFERATIVE:
			file<<"5.0"<<std::endl;
			break;
		case Cell::QUIESCENT:
			file<<"2.0"<<std::endl;
			break;
		case Cell::NECROTIC:
			file<<"3.0"<<std::endl;
			break;
		case Cell::INVASIVE:
			file<<"4.0"<<std::endl;
			break;
		}
	}
	file << "e" << std::endl;
	file.close();

	std::cout<<"Command to plot with Gnuplot: plot 'state"<<std::to_string(stateIteration)<<".dat' using 1:2:3 pt 7 ps 0.5 palette"<<std::endl;
}
