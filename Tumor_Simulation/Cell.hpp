#ifndef CELL_HPP_
#define CELL_HPP_

#include<vector>

class Cell {
public:
	// Cell type definition
	typedef enum {ECM=0, PROLIFERATIVE, QUIESCENT, NECROTIC, INVASIVE} CellType;

	// Constructors
	Cell(double x=0.0, double y=0.0, double z=0.0, double pECM=0.0, bool isBoundary=false, CellType type=ECM, const std::vector<int> &list={});
	Cell(const Cell &other);	// Copy
	Cell(Cell &&other);			// Move

	// Destructor
	~Cell(){}

	// Assign operator
	Cell& operator=(const Cell &other);
	
	// Methods and functions
	double getX() const {return x;}
	double getY() const {return y;}
	double getZ() const {return z;}
	CellType getType(){return type;}
	void setType(CellType t){type = t;}
	double getPECM(){return pECM_density;}
	void setPECM(double p){pECM_density = (p <= 0.0)? 0.0 : p;}
	bool getIsBoundary(){return isBoundary;}
	void setIsBoundary(bool b){isBoundary = b;}
	std::vector<int>& getNeighbors(){return neighbors;}

private:
	double x, y, z;
	CellType type;
	double pECM_density;
	bool isBoundary;
	std::vector<int> neighbors;	// Adjacency list of the cells graph
};

#endif /* CELL_HPP_ */
