#include "hamiltonian.h"
#include "WaveFunction.h"

Hamiltonian::Hamiltonian(int nParticles_, int nDimensions_, double h_, double h2_, int charge_) :

    nDimensions(nDimensions_),
    nParticles(nParticles_),
    h(h_),
    h2(h2_),
    charge(charge_)
{
}


double Hamiltonian::localEnergy(const mat &r, const double &alpha, const double &beta, WaveFunction *function)
{
    double kinEnergy = kineticEnergy(r, alpha, beta, function);
    double potEnergy = potentialEnergy(r);

    return kinEnergy + potEnergy;
}

double Hamiltonian::kineticEnergy(const mat &r, const double &alpha, const double &beta, WaveFunction *function)
{


    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    double waveFunctionCurrent = function->waveFunction(r, alpha, beta);


    // Kinetic energy

    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = function->waveFunction(rMinus, alpha, beta);
            waveFunctionPlus = function->waveFunction(rPlus, alpha, beta);
            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
    kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

    return kineticEnergy;
}


double Hamiltonian::potentialEnergy(const mat &r)
{
    // Potential energy
    double potentialEnergy = 0;
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j)*r(i,j);
        }
        potentialEnergy -= charge / sqrt(rSingleParticle);
    }
    // Contribution from electron-electron potential
//    double r12 = 0;
//    for(int i = 0; i < nParticles; i++) {
//        for(int j = i + 1; j < nParticles; j++) {
//            r12 = 0;
//            for(int k = 0; k < nDimensions; k++) {
//                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
//            }
//            potentialEnergy += 1 / sqrt(r12);
//        }
//    }

    return potentialEnergy;
}

double Hamiltonian::analyticLocalEnergy(const mat &r, const double &alpha, const double &beta)
{
    double r1 = norm(r.row(0),2);
    double r2 = norm(r.row(1),2);
    double r12 = norm((r.row(1) - r.row(0)),2);
    double dot_r12;
    for (int i = 0; i < nDimensions; i++) { dot_r12 += r(0,i)*r(1,i); }
    double energy_l1 = (alpha-charge)*(1/r1 + 1/r2) + 1/r12 - pow(alpha,2);
    double energy_part = ((alpha*(r1+r2))/r12)*(1-dot_r12/(r1*r2)) - 1/(2*pow((1+beta*r12),2)) - 2/r12 + (2*beta)/(1+beta*r12);
    double energy_l2 = energy_l1 + 1/(2*pow((1+beta*r12),2))*energy_part;

    return energy_l2;
}
