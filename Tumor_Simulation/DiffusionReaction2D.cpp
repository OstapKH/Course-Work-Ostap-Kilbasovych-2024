#include <iostream>
#include <fstream>
#include "DiffusionReaction2D.hpp"
#include <boost/numeric/odeint.hpp>

void DiffusionReaction2D::setDefaultDrugsToBoundaryCells() {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            domain[i][j] = (i == 0 || i == height - 1 || j == 0 || j == width - 1) ? diffParam.fi0 : 0;
        }
    }
}

std::vector<std::vector<double>> &DiffusionReaction2D::computeNextState() {

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (i == 0 || i == height - 1 || j == 0 || j == width - 1) {
                double fi0 = calculateFi0(diffParam.Dt * stateIteration / (3600 * 24), 43200);
                domainNext[i][j] = fi0;
            }
        }
    }

    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            double currentValue = domain[i][j];
            double laplacianValue = (laplacian(domain, i, j));

            domainNext[i][j] = currentValue + diffParam.Dt * (
                    (diffParam.D0 * laplacianValue - diffParam.Kmet * currentValue -
                     diffParam.lambda14 * diffParam.nCA *
                     (currentValue / (currentValue + calculateFi0(diffParam.Dt * stateIteration / (3600 * 24), 43200))))
            );
        }
    }
    stateIteration++;
    domain = domainNext;
    return domainNext;
}

void DiffusionReaction2D::storeLastStateToFile() {
    std::ofstream file("state_diffusion_reaction" + std::to_string(stateIteration) + ".dat");
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            double value = domainNext[i][j];
            file << i << " " << j << " " << value << std::endl;
        }
    }
    file << "e" << std::endl;
    file.close();
}

void DiffusionReaction2D::printLastStateToConsole() {
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            double value = domain[i][j];
            std::cout << "\t\t" << value << "\t\t";
        }
        std::cout << std::endl;
    }
}


double DiffusionReaction2D::laplacian(std::vector<std::vector<double>> &arr, int i, int j) {
    double i1_j = arr[i + 1][j];
    double im1_j = arr[i - 1][j];
    double i_j1 = arr[i][j + 1];
    double i_jm1 = arr[i][j - 1];
    double i_j = arr[i][j];
    double Dx2 = diffParam.Dx * diffParam.Dx;
    double result = (i1_j + im1_j + i_j1 + i_jm1 - 4 * i_j) / (Dx2);
    return result;
}

void DiffusionReaction2D::initialiseDomains() {
    domain = std::vector<std::vector<double>>(height, std::vector<double>(width, 0.0));
    domainNext = std::vector<std::vector<double>>(height, std::vector<double>(width, 0.0));
}

// PBPK model
typedef std::vector<double> state_type;

const double k10 = 13.27, k12 = 0.276, k21 = 1.48, V1 = 4.84E3, V2 = 8.0E3; // PBPK model parameters

void odt(const state_type &x, state_type &dx, double t) {
    dx[0] = ((k21 * x[1] * V2) / V1) - ((k12 + k10) * x[0]) + (0.0 / V1);
    dx[1] = ((k12 * x[0] * V1) / V2) - (k21 * x[1]);
}

double DiffusionReaction2D::calculateFi0(double t, double dosagePeriod) {
    const double dt = 0.01;
    boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper;
    state_type x(2);
    x[0] = 10000;
    x[1] = 0.1;

    double yt = 0.0;

    double tWithDosagePeriod = fmod(t, dosagePeriod);

    while (yt < tWithDosagePeriod) {
        stepper.do_step(odt, x, t, dt);
        yt += dt;
    }

    if (this->diffParam.vascular) {
        return log(x[0]) / 10;
    } else {
        return log(x[1]) / 10;
    }
}


void DiffusionReaction2D::fillCellsDiffusionGridWithCells(std::vector<Cell> tumorCells) {
    for (int i = 0; i < tumorCells.size(); i++) {
        cellsDiffusionGrid[floor(tumorCells[i].getX())][floor(tumorCells[i].getY())].push_back(i);
    }
}

void DiffusionReaction2D::calculateAverageVectors(std::vector<Cell> tumorCells) {
    double tempSum = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int c: cellsDiffusionGrid[i][j]) {
                tempSum += tumorCells[c].getPECM();
            }
            averageRo_ECMValues[i][j] = tempSum / (double) cellsDiffusionGrid[i][j].size();
            tempSum = 0;
        }
    }
}


