#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>

using namespace arma;
using namespace std;


class WaveFunction {

public:

    WaveFunction(int nP, int nD);
    double WaveFunction::waveFunction(const mat &r);

private:

    double WaveFunction::jastrowFactor(const mat &r);

    int nDimensions;
    int nParticles;


};

#endif // WAVEFUNCTION_H
