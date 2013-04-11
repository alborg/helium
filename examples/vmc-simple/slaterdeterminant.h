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
    double getDeterminant();
    void updateDeterminant(const mat &r, int i, double &alpha_, double &beta_);
    double beryllium(const mat &r, double &alpha);

private:

    WaveFunction *function;
    int nDimensions;
    int nParticles;
    double alpha;
    mat slaterMatrixUp;
    mat slaterMatrixDown;
    mat slaterMatrixUpNew;
    mat slaterMatrixDownNew;
    mat invSlaterMatrixUp;
    mat invSlaterMatrixDown;

};

#endif // SLATERDETERMINANT_H
