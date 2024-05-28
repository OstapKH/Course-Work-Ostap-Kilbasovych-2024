#include <iostream>
#include <fstream>
#include "DiffusionReaction3D.hpp"
#include <boost/numeric/odeint.hpp>
#include <mpi.h>

DiffusionReaction3D::DiffusionReaction3D(int rank, int size, int domainSize, int numIterations, bool isVascular) {
    this->height = this->width = this->depth = domainSize;
    this->diffParam.vascular = isVascular;
    this->stateIteration = 0;

    initialiseDomains();
    setDefaultDrugsToBoundaryCells();

    MPI_Bcast(&numIterations, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(&domain[0], domainSize*domainSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Calculate part of the domain
        computeNextStateMPI(rank, size);

        std::vector<std::vector<std::vector<double>>> receivedDomain(domain.size(),
                                                                     std::vector<std::vector<double>>(domain[0].size(),
                                                                                                      std::vector<double>(
                                                                                                              domain[0][0].size(),
                                                                                                              0)));
        for (int i = 1; i < size; i++) {
            // Vector to store the part of received newDomain
            int subHeight = int(height) / size;
            MPI_Recv(&receivedDomain[0], domain.size() * domain[0].size(), MPI_DOUBLE, i, 0,
                     MPI_COMM_WORLD, nullptr);
            for (int j = rank * subHeight; j < ((rank + 1) * subHeight); ++j) {
                for (int k = rank * subHeight; k < ((rank + 1) * subHeight); ++k) {
                    for (int l = rank * subHeight; l < ((rank + 1) * subHeight); ++l) {
                        domainNext[j][k][l] = receivedDomain[j][k][l];
                    }
                }
            }
        }

        domain = domainNext;
        storeLastStateToFile();
    } else {
        computeNextStateMPI(rank, size);
        MPI_Send(&domainNext[0], domain.size() * domain[0].size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Bcast(&domain[0], domain.size() * domain[0].size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}


void DiffusionReaction3D::setDefaultDrugsToBoundaryCells() {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < depth; k++) {
                domain[i][j][k] = (i == 0 || i == height - 1 || j == 0 || j == width - 1 || k == 0 || k == depth - 1)
                                  ? diffParam.fi0 : 0;
            }
        }
    }
}

std::vector<std::vector<std::vector<double>>> &DiffusionReaction3D::computeNextState() {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < depth; k++) {
                double fi0 = calculateFi0(diffParam.Dt * stateIteration / (3600 * 24), 43200);
                if (i == 0 || i == height - 1 || j == 0 || j == width - 1 || k == 0 || k == depth - 1) {
                    domainNext[i][j][k] = fi0;
                }
            }
        }
    }

    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            for (int k = 1; k < depth - 1; k++) {
                double currentValue = domain[i][j][k];
                double laplacianValue = (laplacian(domain, i, j, k));
                domainNext[i][j][k] = currentValue + diffParam.Dt * (
                        (diffParam.D0 * laplacianValue - diffParam.Kmet * currentValue -
                         diffParam.lambda14 * diffParam.nCA *
                         (currentValue /
                          (currentValue + calculateFi0(diffParam.Dt * stateIteration / (3600 * 24), 43200))))
                );
            }
        }
    }
    stateIteration++;
    domain = domainNext;
    return domainNext;
}

void DiffusionReaction3D::setBorderValue() {
    std::cout << "Setting border values..." << std::endl;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < depth; k++) {
                double fi0 = calculateFi0(diffParam.Dt * stateIteration / (3600 * 24), 43200);
                if (i == 0 || i == height - 1 || j == 0 || j == width - 1 || k == 0 || k == depth - 1) {
                    domainNext[i][j][k] = fi0;
                }
            }
        }
    }
    std::cout << "Border values set" << std::endl;
}

void DiffusionReaction3D::computeNextStateMPI(int rank, int size) {
    int subHeight = int(height) / size;
    int startHeight = rank * subHeight;
    int endHeight = ((rank + 1) * subHeight);

    if (rank == size - 1) {
        endHeight = int(height); // Last process takes the remainder
    }

    if(rank == 0){
        for (int i = startHeight+1; i < endHeight - 1; ++i) {
            for (int j = startHeight+1; j < endHeight - 1; ++j) {
                for (int k = startHeight+1; k < endHeight - 1; ++k) {
                    double currentValue = domain[i][j][k];
                    double laplacianValue = (laplacian(domain, i, j, k));
                    domainNext[i][j][k] = currentValue + diffParam.Dt * (
                            (diffParam.D0 * laplacianValue - diffParam.Kmet * currentValue -
                             diffParam.lambda14 * diffParam.nCA *
                             (currentValue /
                              (currentValue + calculateFi0(diffParam.Dt * stateIteration / (3600 * 24), 43200))))
                    );
                }
            }
        }
    }
    else{
        for (int i = startHeight; i < endHeight - 1; ++i) {
            for (int j = startHeight; j < endHeight - 1; ++j) {
                for (int k = startHeight; k < endHeight - 1; ++k) {
                    double currentValue = domain[i][j][k];
                    double laplacianValue = (laplacian(domain, i, j, k));
                    domainNext[i][j][k] = currentValue + diffParam.Dt * (
                            (diffParam.D0 * laplacianValue - diffParam.Kmet * currentValue -
                             diffParam.lambda14 * diffParam.nCA *
                             (currentValue /
                              (currentValue + calculateFi0(diffParam.Dt * stateIteration / (3600 * 24), 43200))))
                    );
                }
            }
        }
    }



    setBorderValue();
}

void DiffusionReaction3D::storeLastStateToFile() {
    std::ofstream file("state_diffusion_reaction" + std::to_string(stateIteration) + ".dat");
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            for (int k = 0; k < depth; ++k) {
                double value = domainNext[i][j][k];
                file << i << " " << j << " " << k << " " << value << std::endl;
            }
        }
    }
    file << "e" << std::endl;
    file.close();
}

void DiffusionReaction3D::printLastStateToConsole() {
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            for (int k = 0; k < depth; ++k) {
                double value = domain[i][j][k];
                std::cout << "\t\t" << value << "\t\t";
            }
            std::cout << std::endl;
        }
    }
}

double DiffusionReaction3D::laplacian(std::vector<std::vector<std::vector<double>>> &arr, int i, int j, int k) {
    double i1_j_k = arr[i + 1][j][k];
    double im1_j_k = arr[i - 1][j][k];
    double i_j1_k = arr[i][j + 1][k];
    double i_jm1_k = arr[i][j - 1][k];
    double i_j_k1 = arr[i][j][k + 1];
    double i_j_km1 = arr[i][j][k - 1];
    double i_j_k = arr[i][j][k];
    double Dx2 = diffParam.Dx * diffParam.Dx;
    double result = (i1_j_k + im1_j_k + i_j1_k + i_jm1_k + i_j_k1 + i_j_km1 - 6 * i_j_k) / (Dx2);
    return result;
}

void DiffusionReaction3D::initialiseDomains() {
    domain = std::vector<std::vector<std::vector<double>>>(
            height, std::vector<std::vector<double>>(width, std::vector<double>(depth, 0.0)));
    domainNext = std::vector<std::vector<std::vector<double>>>(
            height, std::vector<std::vector<double>>(width, std::vector<double>(depth, 0.0)));
}

// PBPK model

typedef std::vector<double> state_type;

const double k10 = 13.27, k12 = 0.276, k21 = 1.48, V1 = 4.84E3, V2 = 8.0E3; // PBPK model parameters

void odt3D(const state_type &x, state_type &dx, double t) {
    dx[0] = ((k21 * x[1] * V2) / V1) - ((k12 + k10) * x[0]) + (0.0 / V1);
    dx[1] = ((k12 * x[0] * V1) / V2) - (k21 * x[1]);
}

double DiffusionReaction3D::calculateFi0(double t, double dosagePeriod) {
    const double dt = 0.01;
    boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper;
    state_type x(2);
    x[0] = 10000;
    x[1] = 0.1;

    double yt = 0.0;

    double tWithDosagePeriod = fmod(t, dosagePeriod);

    while (yt < tWithDosagePeriod) {
        stepper.do_step(odt3D, x, t, dt);
        yt += dt;
    }

    if (this->diffParam.vascular) {
        return log(x[0]) / 10;
    } else {
        return log(x[1]) / 10;
    }
}

void DiffusionReaction3D::fillCellsDiffusionGridWithCells(std::vector<Cell> tumorCells) {
    for (int i = 0; i < tumorCells.size(); i++) {
        cellsDiffusionGrid[floor(tumorCells[i].getX())][floor(tumorCells[i].getY())][floor(
                tumorCells[i].getZ())].push_back(i);
    }
}

void DiffusionReaction3D::calculateAverageVectors(std::vector<Cell> tumorCells) {
    double tempSum = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < depth; k++) {
                for (int c: cellsDiffusionGrid[i][j][k]) {
                    tempSum += tumorCells[c].getPECM();
                }
                averageRo_ECMValues[i][j][k] = tempSum / (double) cellsDiffusionGrid[i][j][k].size();
                tempSum = 0;
            }
        }
    }
}

std::vector<double> flatten3D(const std::vector<std::vector<std::vector<double>>> &vec) {
    std::vector<double> flat;
    for (const auto &matrix : vec) {
        for (const auto &row : matrix) {
            flat.insert(flat.end(), row.begin(), row.end());
        }
    }
    return flat;
}

void unflatten3D(const std::vector<double> &flat, std::vector<std::vector<std::vector<double>>> &vec, int height, int width, int depth) {
    int index = 0;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            for (int k = 0; k < depth; ++k) {
                vec[i][j][k] = flat[index++];
            }
        }
    }
}

