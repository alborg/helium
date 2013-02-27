#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "WaveFunction.h"


class Hamiltonian
{
public:
    Hamiltonian(int nParticles_, int nDimensions_, double h_, double h2_, int charge_);
    double localEnergy(const mat &r, const double &alpha, const double &beta, WaveFunction *function);
    double analyticLocalEnergy(const mat &r, const double &alpha, const double &beta);

private:
    double kineticEnergy(const mat &r, const double &alpha, const double &beta, WaveFunction *function);
    double potentialEnergy(const mat &r);


    int nDimensions;
    int nParticles;
    double h;
    double h2;
    int charge;

};

#endif // HAMILTONIAN_H
