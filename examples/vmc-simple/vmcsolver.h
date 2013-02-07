#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>

using namespace arma;

class VMCSolver
{
public:
    VMCSolver();

    mat runMonteCarloIntegration(const double &alpha_in, const double &beta_in);

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

    int nCycles;

    double alpha;
    double beta;

    mat rOld;
    mat rNew;
};

#endif // VMCSOLVER_H
