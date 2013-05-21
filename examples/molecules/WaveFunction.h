#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>

using namespace arma;
using namespace std;


class WaveFunction {

public:


    WaveFunction(int &nDimensions_, int &nProtons_, int &nElectrons_);

    double waveFunction(const mat &r, const mat &rProtons, double alpha, double beta);

private:

    double jastrowFactor(const mat &r, double beta);

    int nDimensions;
    int nProtons;
    int nElectrons;
    int nParticles;
    double alpha;
    double beta;

};

#endif // WAVEFUNCTION_H
