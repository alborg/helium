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
    vec gradientWaveFunction(const mat &r, int i, double ratio, double alpha, double beta);
    double laPlaceWaveFunction(const mat &r, double alpha, double beta);
    double laPlaceWaveFunctionNum(const mat &r, double alpha, double beta);
    double getDeterminant();
    double getInvDeterminant();
    vec getStates(const mat &r, int i, double rtot, double alpha, double beta);
    void updateDeterminant(const mat &rNew, const mat &rOld, int i, double &alpha_, double &beta_, double ratio);
    double beryllium(const mat &r, double &alpha);
    vec gradientWaveFunctionNum(const mat &r, int i, double alpha_, double beta_);
    double getRatioDeterminantNum(int i, const mat &rOld, const mat &rNew, double alpha, double beta);

private:

    WaveFunction *function;
    int nDimensions;
    int nParticles;
    mat slaterMatrixUp;
    mat slaterMatrixDown;
    mat invSlaterMatrixUp;
    mat invSlaterMatrixDown;

};

#endif // SLATERDETERMINANT_H
