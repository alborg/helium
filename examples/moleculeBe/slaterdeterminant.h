#ifndef SLATERDETERMINANT_H
#define SLATERDETERMINANT_H

#include <armadillo>
#include "WaveFunction.h"


using namespace arma;
using namespace std;

class slaterDeterminant
{
public:
    slaterDeterminant(int nDimensions_, int nProtons_, int nElectrons_);
    void buildDeterminant(const mat &r, double &alpha_);
    double getRatioDeterminant(int i, const mat &r, double alpha);
    vec gradientWaveFunction(const mat &r, int i, double ratio, double alpha);
    double laPlaceWaveFunction(const mat &r, int i, double alpha);
    double laPlaceWaveFunctionNum(const mat &r, double alpha);
    double getDeterminant();
    double getInvDeterminant();
    vec getStates(double rtot, double alpha);
    void updateDeterminant(const mat &rNew, const mat &rOld, int i, double &alpha_, double ratio);
    double beryllium(const mat &r, double &alpha);
    vec gradientWaveFunctionNum(const mat &r, int i, double alpha_);
    double getRatioDeterminantNum(int i, const mat &rOld, const mat &rNew, double alpha);
    double dWaveFunction_dalpha(const mat &r, double alpha);

private:

    WaveFunction *function;
    int nDimensions;
    int nElectrons;
    int nProtons;
    int nParticles;
    mat slaterMatrixUp1;
    mat slaterMatrixDown1;
    mat slaterMatrixUp2;
    mat slaterMatrixDown2;
    mat invSlaterMatrixUp1;
    mat invSlaterMatrixDown1;
    mat invSlaterMatrixUp2;
    mat invSlaterMatrixDown2;

};

#endif // SLATERDETERMINANT_H
