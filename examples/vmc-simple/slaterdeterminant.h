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
    void buildDeterminant(const mat &r, double &alpha_);
    double getDeterminant(const mat &r, double &alpha_);
    double beryllium(const mat &r, double &alpha);

private:

    int nDimensions;
    int nParticles;
    double alpha;
    mat slaterMatrixUp;
    mat slaterMatrixDown;
    WaveFunction *function;

};

#endif // SLATERDETERMINANT_H
