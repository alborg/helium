#ifndef SLATERDETERMINANT_H
#define SLATERDETERMINANT_H

#include <armadillo>


using namespace arma;
using namespace std;

class slaterDeterminant
{
public:
    slaterDeterminant(int nParticles_, int nDimensions_);
    double beryllium(const mat &r, double &alpha);

private:
    double psi1s(double r);
    double psi2s(double r);

    int nParticles;
    int nDimensions;
    double alpha;
};

#endif // SLATERDETERMINANT_H
