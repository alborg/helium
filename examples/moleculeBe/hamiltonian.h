#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "slaterdeterminant.h"
#include "correlation.h"


class Hamiltonian
{
public:
    Hamiltonian(int nDimensions_, double h_, double h2_, int nProtons_, int nElectrons_, double R_, const mat &rProtons_);
    double localEnergy(const mat &r, const double &alpha, const double &beta, slaterDeterminant *slater, correlation *corr);
    vec dPsi(const mat &r, double alpha, double beta, slaterDeterminant *slater, correlation *corr);
    double analyticEnergyHe(const mat &r, const double &alpha, const double &beta);
    double analyticdEnergyHe(const mat &r, const double &alpha, const double &beta);

private:
    double kineticEnergy(const mat &r, const double &alpha, const double &beta, slaterDeterminant *slater, correlation *corr);
    double potentialEnergy(const mat &r);


    int nDimensions;
    int nProtons;
    int nElectrons;
    int nParticles;
    double h;
    double h2;
    double R;
    mat rProtons;

};

#endif // HAMILTONIAN_H
