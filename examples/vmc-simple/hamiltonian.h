#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "WaveFunction.h"


class Hamiltonian
{
public:
    Hamiltonian(int nParticles_, int nDimensions_, double h_, double h2_, int charge_);
    double localEnergy(const mat &r, const double &alpha, const double &beta, WaveFunction *function);

private:
    int nDimensions;
    int nParticles;
    int h;
    int h2;
    double charge;

};

#endif // HAMILTONIAN_H
