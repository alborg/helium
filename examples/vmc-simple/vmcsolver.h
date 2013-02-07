#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>

using namespace arma;

class VMCSolver
{
public:
    VMCSolver();

    mat runMonteCarloIntegration(const double &alpha, const double &beta);

private:
    double waveFunction(const mat &r);
    double localEnergy(const mat &r);

    int nDimensions;
    int charge;
    double stepLength;
    int nParticles;

    double h;
    double h2;

    long idum;

    double alpha;
    double beta;

    int accepted_steps;

    int nCycles;

    mat rOld;
    mat rNew;
};

#endif // VMCSOLVER_H
