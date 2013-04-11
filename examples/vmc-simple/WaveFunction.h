#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
#include "slaterdeterminant.h"

using namespace arma;
using namespace std;


class WaveFunction {

public:

    WaveFunction(int &nParticles_, int &nDimensions_, slaterDeterminant *slater_);

    double waveFunction(const mat &r, double alpha, double beta);
    double gradientWaveFunction(const mat &r, int particle, int dimension, double alpha, double beta);
    double psi1s(double r, double alpha);
    double psi2s(double r, double alpha);
    double psi2p(double r, double alpha);
    double dPsi1s(double rtot, double variable, double alpha);
    double dPsi2s(double rtot, double variable, double alpha);
    double dPsi2p(double rtot, double variable, double alpha);

private:

    double jastrowFactor(const mat &r, double beta);

    int nDimensions;
    int nParticles;
    double alpha;
    double beta;
    slaterDeterminant *slater;

};

#endif // WAVEFUNCTION_H
