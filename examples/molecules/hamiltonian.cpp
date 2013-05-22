#include "hamiltonian.h"
#include "WaveFunction.h"

Hamiltonian::Hamiltonian(int nProtons_, int nElectrons_, int nDimensions_, double h_, double h2_, int charge_) :

    nDimensions(nDimensions_),
    nProtons(nProtons_),
    nElectrons(nElectrons_),
    nParticles(nProtons*nElectrons),
    h(h_),
    h2(h2_),
    charge(charge_)
{
}



//Find the local energy (expectation value of the energy) numerically
double Hamiltonian::localEnergy(double R, const mat &r, const mat &rProtons, const double &alpha, const double &beta, WaveFunction *function)
{


    double kinEnergy = kineticEnergy(r, rProtons, alpha, beta, function);
    double potEnergy = potentialEnergy(R, r, rProtons);
    //cout<<"kinetic, potential: "<<kinEnergy<<" "<<potEnergy<<endl;

    return kinEnergy + potEnergy;
}



//Find the kinetic energy part of the local energy
double Hamiltonian::kineticEnergy(const mat &r, const mat rProtons, const double &alpha, const double &beta, WaveFunction *function)
{

    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    double waveFunctionCurrent = function->waveFunction(r, rProtons, alpha, beta); //Find wavefunction for r

    //Second derivative (del^2):

    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = function->waveFunction(rMinus, rProtons, alpha, beta);
            waveFunctionPlus = function->waveFunction(rPlus, rProtons, alpha, beta);
            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
    kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

    return kineticEnergy;

}


double Hamiltonian::potentialEnergy(double R, const mat &r, const mat &rProtons)
{

    double potentialE = 0;

    //Contribution from electron - proton potential (1/rep)
    double rp = 0;

    for(int e=0; e<nParticles; e++) {
        for(int p=0; p<nProtons; p++) {
            rp = 0;
            for(int d=0; d<nDimensions; d++) rp += (r(e,d) - rProtons(p,d))*(r(e,d) - rProtons(p,d));
            potentialE -= charge/sqrt(rp);
        }
    }



    // Contribution from electron-electron potential (1/rij part)
    double r12 = 0;
    for(int i = 1; i < nParticles; i++) {
        for(int j = 0; j < i; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) r12 += pow((r(i,k) - r(j,k)),2);
            potentialE += 1 / sqrt(r12);
        }
    }

    //Contribution from proton-proton potential 1/R

    potentialE += abs(1/R);



    return potentialE;
}
