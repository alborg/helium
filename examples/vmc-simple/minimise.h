#ifndef MINIMISE_H
#define MINIMISE_H

#include <armadillo>
#include "slaterdeterminant.h"
#include "correlation.h"
#include "hamiltonian.h"

using namespace arma;
using namespace std;

class minimise
{
public:
    minimise(int nP);
    void dfpmin(int var,const mat &r,double var2, slaterDeterminant *slater,correlation *corr,
                vec &p, int n, double gtol, int *iter, double *fret, Hamiltonian *hamiltonian);
    double steepestDescent(vec(*func)(double alpha, double beta, int mpi_steps,
                                        slaterDeterminant *slater, Hamiltonian *hamiltonian, correlation *corr,
                                        double *allEnergies),double alpha, double beta, int min_steps,
                                        slaterDeterminant *slater, Hamiltonian *hamiltonian, correlation *corr);

private:
    void lnsrch(int var,const mat &r,double var2, slaterDeterminant *slater,correlation *corr,
                int n, vec &xold, double fold, vec &g, vec &p, vec &x,
                double *f, double stpmax, int *check, Hamiltonian *hamiltonian);

    int nParticles;
};

#endif // MINIMISE_H
