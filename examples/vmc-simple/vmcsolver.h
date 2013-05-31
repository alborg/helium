#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include "WaveFunction.h"
#include "hamiltonian.h"
#include "slaterdeterminant.h"
#include "correlation.h"
#include "minimise.h"

using namespace arma;
using namespace std;

class VMCSolver
{
public:
    VMCSolver();
    void runMonteCarloIntegration(int argc, char* argv[]);


private:

    vec MCSampling(long idum, double alpha, double beta, int mpi_steps,
                    slaterDeterminant *slater, Hamiltonian *hamiltonian, correlation *corr,
                    double *allEnergies);
    mat quantumForce(const mat &r, double alpha_, double beta_, double wf, WaveFunction *function);
    double gaussianDeviate(long *idum);
    vec MCImportance(long idum, double alpha, double beta, int mpi_steps,
                      slaterDeterminant *slater, Hamiltonian *hamiltonian, correlation *corr,
                      double *allEnergies);
    vec gradE(vec dPsi, double Elocal, vec dPsi_Elocal);
    vec steepest_descent(long idum, vec &pnew, int n, double gtol, int min_steps, int *iter, double *fret,
                           slaterDeterminant *slater, correlation *corr, Hamiltonian *hamiltonian);


    mat rOld;
    mat rNew;

    int nDimensions;
    double h;
    double h2;
    double timestep;
    double D;
    double stepLength;
    int nCycles;
    int charge;
    int nParticles;
    double alpha;
    double beta;

    bool minimise_var;
    int min_steps;

};

#endif // VMCSOLVER_H
