#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <armadillo>
#include "WaveFunction.h"

using namespace arma;

class Hamiltonian
{
public:
    Hamiltonian(int nProtons_, int nElectrons_, int nDimensions_, double h_, double h2_);
    double localEnergy(double R, const mat &r, const mat &rProtons, const double &alpha, const double &beta, WaveFunction *function);

private:
    double kineticEnergy(const mat &r, const mat rProtons, const double &alpha, const double &beta, WaveFunction *function);
    double potentialEnergy(double R, const mat &r, const mat &rProtons);


    int nDimensions;
    int nProtons;
    int nElectrons;
    int nParticles;
    double h;
    double h2;
    int charge;

};

#endif // HAMILTONIAN_H

