#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>

using namespace arma;
using namespace std;


class WaveFunction {

public:

    WaveFunction(int &nParticles_, int &nDimensions_);

    double waveFunction(const mat &r, double alpha, double beta);
    double gradientWaveFunction(const mat &r, int particle, int dimension, double alpha, double beta);
    double psi1s(double r, double alpha);
    double psi2s(double r, double alpha);
    double psi2p(double r, double alpha);
    vec dPsi1s(double rtot, const mat &r, double alpha);
    vec dPsi2s(double rtot, const mat &r, double alpha);
    vec dPsi2p0(double rtot, const mat &r, double alpha);
    vec dPsi2p_1(double rtot, const mat &r, double alpha);
    vec dPsi2p1(double rtot, const mat &r, double alpha);
    double d2Psi1s(double rtot, const mat &r, double alpha);
    double d2Psi2s(double rtot, const mat &r, double alpha);
    double d2Psi2p0(double rtot, const mat &r, double alpha);
    double d2Psi2p_1(double rtot, const mat &r, double alpha);
    double d2Psi2p1(double rtot, const mat &r, double alpha);


private:

    double jastrowFactor(const mat &r, double beta);

    int nDimensions;
    int nParticles;
    double alpha;
    double beta;

};

#endif // WAVEFUNCTION_H
