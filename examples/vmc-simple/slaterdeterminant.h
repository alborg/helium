#ifndef SLATERDETERMINANT_H
#define SLATERDETERMINANT_H

#include <armadillo>


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
    double psi1s(double r);
    double psi2s(double r);
    double psi2p(double r);


    int nDimensions;
    int nParticles;
    double alpha;
    mat slaterMatrixUp;
    mat slaterMatrixDown;

};

#endif // SLATERDETERMINANT_H
