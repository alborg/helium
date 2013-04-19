#ifndef SLATERDETERMINANT_H
#define SLATERDETERMINANT_H

#include <armadillo>
#include "WaveFunction.h"


using namespace arma;
using namespace std;

class slaterDeterminant
{
public:
    slaterDeterminant(int nParticles_, int nDimensions_);
    void buildDeterminant(const mat &r, double &alpha_, double &beta_);
    double getRatioDeterminant(int i, const mat &r, double alpha, double beta);
    double laPlaceWaveFunction(const mat &r, double alpha, double beta);
    double getDeterminant();
    vec getStates(const mat &r, int i, double rtot, double alpha, double beta);
    void updateDeterminant(const mat &rNew, const mat &rOld, int i, double &alpha_, double &beta_, double ratio);
    double beryllium(const mat &r, double &alpha);

private:

    WaveFunction *function;
    int nDimensions;
    int nParticles;
    double alpha;
    mat slaterMatrixUp;
    mat slaterMatrixDown;
    mat invSlaterMatrixUp;
    mat invSlaterMatrixDown;

};

#endif // SLATERDETERMINANT_H
