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
    vec MCImportance(long idum, double alpha, double beta, int mpi_steps,
                      slaterDeterminant *slater, Hamiltonian *hamiltonian, correlation *corr,
                      double *allEnergies);
    vec MCSampling(long idum, double alpha, double beta, int mpi_steps,
                    slaterDeterminant *slater, Hamiltonian *hamiltonian, correlation *corr,
                    double *allEnergies);
    mat quantumForce(const mat &r, double alpha_, double beta_, double wf, WaveFunction *function);
    double gaussianDeviate(long *idum);
    double gradE(vec dPsi, double Elocal, vec dPsi_Elocal, vec &g);
    void dfpmin(vec &p, int n, double gtol, int min_steps, int *iter, double *fret,
                           slaterDeterminant *slater,correlation *corr, Hamiltonian *hamiltonian);
    void lnsrch(long idum, double alpha_start, double beta_start, int n, vec &xold, double fold, vec &g, vec &p, vec &x, double *f, double stpmax,
                           int *check, int min_steps, slaterDeterminant *slater, Hamiltonian *hamiltonian,
                           correlation *corr, double *allEnergies);

    void printFile(const char &file_energies, const char &file_energySquareds, const char &file_alpha, const char &file_sigma, const mat &energies, const mat &energiesSquared, const vec alphas, const vec betas);

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
