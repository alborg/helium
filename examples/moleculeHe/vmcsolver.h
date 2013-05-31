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
    vec MCImportance(long idum, double alpha, double beta, int mpi_steps, Hamiltonian *hamiltonian, WaveFunction *function,
                     double *allEnergies);
    mat quantumForce(const mat &r, const mat &rProtons, double alpha_, double beta_, double wf,
                     WaveFunction *function);
    double gaussianDeviate(long *idum);
    vec steepest_descent(long idum, vec &p, int n, double gtol, int min_steps, int *iter, double *fret,
                           Hamiltonian *hamiltonian, WaveFunction *function);
    vec gradE(vec dPsi, double Elocal, vec dPsi_Elocal);


    mat rOld;
    mat rNew;

    int nDimensions;
    double h;
    double h2;
    double timestep;
    double D;
    double stepLength;
    int nCycles;
    double alpha;
    double beta;
    double R;
    int nProtons;
    int nElectrons;
    int nParticles;
    int charge;

    bool minimise_var;
    int min_steps;
    bool printToFile;

    mat rProtons;


};

#endif // VMCSOLVER_H
