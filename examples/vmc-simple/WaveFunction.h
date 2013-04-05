#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
#include "slaterdeterminant.h"

using namespace arma;
using namespace std;


class WaveFunction {

public:

    WaveFunction(int &nParticles_, int &nDimensions_);
    double waveFunction(const mat &r, double alpha, double beta);
    void buildDeterminant(const mat &r, double alpha);
    double psi1s(double r, double alpha);
    double psi2s(double r, double alpha);
    double psi2p(double r, double alpha);

private:

    double jastrowFactor(const mat &r);

    int nDimensions;
    int nParticles;
    double alpha;
    double beta;
    slaterDeterminant *slater;

};

#endif // WAVEFUNCTION_H
