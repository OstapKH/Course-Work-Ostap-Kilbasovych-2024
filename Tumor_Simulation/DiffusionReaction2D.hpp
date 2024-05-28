#pragma once

#include "Cell.hpp"
#include <vector>

class DiffusionReaction2D {
public:
    // Diffusion-Reaction model simulation parameters type definition
    struct DiffusionReactionParameters {
        double fi0, D0, lambda14, Dt, Dx, Kmet, nCA, dosagePeriod;
        bool vascular;

        DiffusionReactionParameters(double fi0 = 1, double D0 = 1.3e-6, double lambda14 = 2.5e-4, double Dt = 100,
                                    double Dx = 0.1, double Kmet = 2.0e-4, double nCA = 0.15, bool vascular = true, double dosagePeriod = 7200)
                : fi0(fi0), D0(D0), lambda14(lambda14), Dt(Dt), Dx(Dx), Kmet(Kmet), nCA(nCA), vascular(vascular), dosagePeriod(dosagePeriod) {

        }
    };

    // Constructors
    DiffusionReaction2D() : height(10.0), width(10.0), deep(10.0), stateIteration(0) {
        initialiseDomains();
    }

    DiffusionReaction2D(double height, double width, double deep, bool vascular)
            : height(height), width(width), deep(deep), stateIteration(0) {
        this->diffParam.vascular = vascular;
        initialiseDomains();
    }

    DiffusionReaction2D(double sizeOfGrid) : height(sizeOfGrid), width(sizeOfGrid), deep(sizeOfGrid), stateIteration(0) {
        initialiseDomains();
    }

    DiffusionReaction2D(double sizeOfGrid, bool vascular) : height(sizeOfGrid), width(sizeOfGrid), deep(sizeOfGrid),
                                                                                      stateIteration(0) {
        this->diffParam.vascular = vascular;
        initialiseDomains();
    }

    DiffusionReaction2D(double height, double width, double deep, double fi0, double D0,
                      double lambda14, double Dt, double Dx, double Kmet, double nCA, bool vascular, double dosagePeriod)
            : height(height), width(width), deep(deep), stateIteration(0) {
        this->diffParam.fi0 = fi0;
        this->diffParam.D0 = D0;
        this->diffParam.lambda14 = lambda14;
        this->diffParam.Dt = Dt;
        this->diffParam.Dx = Dx;
        this->diffParam.Kmet = Kmet;
        this->diffParam.nCA = nCA;
        this->diffParam.vascular = vascular;
        this->diffParam.dosagePeriod = dosagePeriod;
        initialiseDomains();
    }

    // Destructor
    ~DiffusionReaction2D() {}

    // Methods and functions
    void setDefaultDrugsToBoundaryCells();

    void storeLastStateToFile();

    std::vector<std::vector<double>> &getLastState() { return domainNext; }

    void calculateAverageRoECMDiffusion(std::vector<Cell> tumorCells);

    std::vector<std::vector<double>> &computeNextState();

    void printLastStateToConsole();

    std::vector<std::vector<double>> &getAvgP_ECM_density() { return averageRo_ECMValues; }


private:
    // Private methods
    double laplacian(std::vector<std::vector<double>> &arr, int i, int j);

    double calculateFi0(double t, double dosagePeriod);

    // Attributes
    double height, width, deep; // Space size
    DiffusionReactionParameters diffParam; // Diffusion-Reaction model simulation parameters
    std::vector<std::vector<double>> domain; // Old state of the model
    std::vector<std::vector<double>> domainNext; // New/last state of the model
    std::vector<std::vector<double>> auxPtr; // Auxiliary vector used to change between old and new states
    int stateIteration; // Counts how many times the computeNextState method has been called
    void initialiseDomains();

    std::vector<std::vector<std::vector<int>>> cellsDiffusionGrid;

    std::vector<std::vector<double>> averageRo_ECMValues;

    void fillCellsDiffusionGridWithCells(std::vector<Cell> tumorCells);

    void calculateAverageVectors(std::vector<Cell> tumorCells);
};
