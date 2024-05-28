#pragma once

#include <vector>
#include "Cell.hpp"
#include <mpi.h>

class DiffusionReaction3D {
public:
    DiffusionReaction3D(int rank, int size, int domainSize, int numIterations, bool isVascular);

    // Diffusion-Reaction model simulation parameters type definition
    struct DiffusionReactionParameters {
        double fi0, D0, lambda14, Dt, Dx, Kmet, nCA;
        bool vascular;

        DiffusionReactionParameters(double fi0 = 1, double D0 = 1.3e-6, double lambda14 = 2.5e-4, double Dt = 100,
                                    double Dx = 1, double Kmet = 2.0e-4, double nCA = 0.15, bool vascular = true)
                : fi0(fi0), D0(D0), lambda14(lambda14), Dt(Dt), Dx(Dx), Kmet(Kmet), nCA(nCA), vascular(vascular) {

        }
    };

    // Constructors
    DiffusionReaction3D() : height(5.0), width(5.0), depth(5.0), stateIteration(0) {
        initialiseDomains();
    }

    DiffusionReaction3D(double height, double width, double depth, bool vascular)
            : height(height), width(width), depth(depth), stateIteration(0) {
        this->diffParam.vascular = vascular;
        initialiseDomains();
    }

    DiffusionReaction3D(double sizeOfGrid) : height(sizeOfGrid), width(sizeOfGrid), depth(sizeOfGrid), stateIteration(0) {
        initialiseDomains();
    }

    DiffusionReaction3D(double sizeOfGrid, bool vascular) : height(sizeOfGrid), width(sizeOfGrid), depth(sizeOfGrid),
                                                          stateIteration(0) {
        this->diffParam.vascular = vascular;
        initialiseDomains();
    }

    DiffusionReaction3D(int height, int width, int depth)
            : height(height), width(width), depth(depth), stateIteration(0) {
        initialiseDomains();
    }


    DiffusionReaction3D(double height, double width, double depth, double fi0, double D0,
                      double lambda14, double Dt, double Dx, double Kmet, double nCA, bool vascular)
            : height(height), width(width), depth(depth), stateIteration(0) {
        this->diffParam.fi0 = fi0;
        this->diffParam.D0 = D0;
        this->diffParam.lambda14 = lambda14;
        this->diffParam.Dt = Dt;
        this->diffParam.Dx = Dx;
        this->diffParam.Kmet = Kmet;
        this->diffParam.nCA = nCA;
        this->diffParam.vascular = vascular;
        initialiseDomains();
    }

    // Destructor
    ~DiffusionReaction3D() {}

    // Methods and functions
    void setDefaultDrugsToBoundaryCells();

    void storeLastStateToFile();

    std::vector<std::vector<std::vector<double>>> &getLastState() { return domainNext; }

    void calculateAverageRoECMDiffusion(std::vector<Cell> tumorCells);

    std::vector<std::vector<std::vector<double>>> &computeNextState();

    void computeNextStateMPI(int rank, int size);


    void printLastStateToConsole();

    std::vector<std::vector<std::vector<double>>> &getAvgP_ECM_density() { return averageRo_ECMValues; }


private:
    // Private methods
    double laplacian(std::vector<std::vector<std::vector<double>>> &arr, int i, int j, int k);

    double calculateFi0(double t, double dosagePeriod);

    // Attributes
    double height, width, depth; // Space size
    DiffusionReactionParameters diffParam;    // Diffusion-Reaction model simulation parameters
    std::vector<std::vector<std::vector<double>>> domain;   // Old state of the model
    std::vector<std::vector<std::vector<double>>> domainNext; // New/last state of the model
    std::vector<std::vector<std::vector<double>>> auxPtr;   // Auxiliar vector used to change between old and new states
    int stateIteration; // Counts how many times the computeNextState method has been called
    void initialiseDomains();

    std::vector<std::vector<std::vector<std::vector<int>>>> cellsDiffusionGrid;

    std::vector<std::vector<std::vector<double>>> averageRo_ECMValues;

    void fillCellsDiffusionGridWithCells(std::vector<Cell> tumorCells);

    void calculateAverageVectors(std::vector<Cell> tumorCells);

    void setBorderValue();

    std::vector<double> flatten3D(const std::vector<std::vector<std::vector<double>>> &vec);

    void
    unflatten3D(const std::vector<double> &flat, std::vector<std::vector<std::vector<double>>> &vec, int height,
                int width,
                int depth);
};
