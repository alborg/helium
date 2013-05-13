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
    vec MCImportance(double alpha, double beta, int mpi_steps,
                      slaterDeterminant *slater, Hamiltonian *hamiltonian, correlation *corr,
                      double *allEnergies);
    vec MCSampling(double alpha, double beta, int mpi_steps,
                    slaterDeterminant *slater, Hamiltonian *hamiltonian, correlation *corr,
                    double *allEnergies);
    mat quantumForce(const mat &r, double alpha_, double beta_, double wf, WaveFunction *function);
    double gaussianDeviate(long *idum);

    void printFile(const char &file_energies, const char &file_energySquareds, const char &file_alpha, const char &file_sigma, const mat &energies, const mat &energiesSquared, const vec alphas, const vec betas);

    mat rOld;
    mat rNew;

    int nDimensions;
    double h;
    double h2;
    long idum;
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
