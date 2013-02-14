#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>

using namespace arma;
using namespace std;


class WaveFunction {

public:

    WaveFunction(int nParticles_, int nDimensions_, double alpha_, double beta_);
    double waveFunction(const mat &r);

private:

    double jastrowFactor(const mat &r);

    int nDimensions;
    int nParticles;
    double alpha;
    double beta;


};

#endif // WAVEFUNCTION_H
