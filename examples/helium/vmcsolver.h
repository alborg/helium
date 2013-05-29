#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include "WaveFunction.h"
#include "hamiltonian.h"

using namespace arma;
using namespace std;

class VMCSolver
{
public:
    VMCSolver();
    void runMonteCarloIntegration(int argc, char* argv[]);

private:
    void MCImportance(double alpha, double beta, int mpi_steps, WaveFunction *function, Hamiltonian *hamiltonian, double &energySum, double &energySquaredSum, double *allEnergies);
    void MCSampling(double alpha, double beta, int mpi_steps, WaveFunction *function, Hamiltonian *hamiltonian, double &energySum, double &energySquaredSum, double *allEnergies);
    mat quantumForce(const mat &r, double alpha_, double beta_, double wf, WaveFunction *function);
    double gaussianDeviate(long *idum);

    void printFile(const char &file_energies, const char &file_energySquareds, const char &file_alpha, const mat &energies, const mat &energiesSquared, const vec alphas, const vec betas);

    int nDimensions;
    int charge;
    int nParticles;

    double h;
    double h2;

    long idum;

    int nCycles;

    mat rOld;
    mat rNew;

    double alpha_min;
    double alpha_max;
    int alpha_steps;
    double beta_min;
    double beta_max;
    int beta_steps;
    double timestep;
    double D;
    double stepLength;


};

#endif // VMCSOLVER_H
