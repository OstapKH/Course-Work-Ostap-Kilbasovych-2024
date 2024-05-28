#include "Cell.hpp"

Cell::Cell(double x, double y, double z, double pECM, bool isBoundary, CellType type, const std::vector<int> &list): x(x), y(y), z(z), type(type), pECM_density(pECM), isBoundary(isBoundary), neighbors(std::vector<int>(list)){}

Cell::Cell(const Cell &other): x(other.x), y(other.y), z(other.z), type(other.type), pECM_density(other.pECM_density), isBoundary(other.isBoundary), neighbors(std::vector<int>(other.neighbors)){}

Cell::Cell(Cell &&other): x(other.x), y(other.y), z(other.z), type(other.type), pECM_density(other.pECM_density), isBoundary(other.isBoundary){
	neighbors = std::move(other.neighbors);
}

Cell& Cell::operator=(const Cell &other){
	if(&other != this){
		x = other.x;
		y = other.y;
		z = other.z;
		type = other.type;
		pECM_density = other.pECM_density;
		isBoundary = other.isBoundary;
		neighbors.clear();
		neighbors = std::vector<int>(other.neighbors);
	}
	return *this;
}

