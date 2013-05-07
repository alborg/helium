#ifndef CORRELATION_H
#define CORRELATION_H
#include <armadillo>

using namespace arma;
using namespace std;

class correlation
{
public:
    correlation(int nParticles_, int nDimensions_);
    double jastrow(const mat &r, double beta);
    double jastrowNum(const mat &r, double beta);
    double getRatioJastrow(int i, const mat &rOld, const mat &rNew, double beta);
    vec gradientWaveFunction(const mat &r, int i, double beta);
    vec gradientWaveFunctionNum(const mat &r, int i, double beta);
    double laPlaceWaveFunction(const mat &r, double beta);
    double laPlaceWaveFunctionNum(const mat &r, double beta);


private:
    int nParticles;
    int nDimensions;
};

#endif // CORRELATION_H
