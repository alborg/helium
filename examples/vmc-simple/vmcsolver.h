#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include "WaveFunction.h"

using namespace arma;
using namespace std;

class VMCSolver
{
public:
    VMCSolver();

    void runMonteCarloIntegration();

private:
    void printFile(const char &file_energies, const char &file_energySquareds, const char &file_alpha, const mat &energies, const mat &energiesSquared, const vec alphas, const vec betas);

    int nDimensions;
    int charge;
    double stepLength;
    int nParticles;

    double h;
    double h2;

    long idum;

    int nCycles;

    mat rOld;
    mat rNew;

    double alpha_min;
    double alpha_max;
    int alpha_steps;
    double beta_min;
    double beta_max;
    int beta_steps;



};

#endif // VMCSOLVER_H
