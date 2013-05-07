#include "hamiltonian.h"
#include "slaterdeterminant.h"
#include "correlation.h"

Hamiltonian::Hamiltonian(int nParticles_, int nDimensions_, double h_, double h2_, int charge_) :

    nDimensions(nDimensions_),
    nParticles(nParticles_),
    h(h_),
    h2(h2_),
    charge(charge_)
{
}


//Find the local energy (expectation value of the energy) numerically
double Hamiltonian::localEnergy(const mat &r, const double &alpha, const double &beta, slaterDeterminant *slater, correlation *corr)
{


    double kinEnergy = kineticEnergy(r, alpha, beta, slater,corr);
    double potEnergy = potentialEnergy(r);

    return kinEnergy + potEnergy;
}

//Find the kinetic energy part of the local energy
double Hamiltonian::kineticEnergy(const mat &r, const double &alpha, const double &beta, slaterDeterminant *slater, correlation *corr)
{

    vec gradSlater = zeros<vec>(nDimensions,1);
    vec gradCorr = zeros<vec>(nDimensions,1);
    double gradProduct = 0;

    double laPlaceSlater = slater->laPlaceWaveFunction(r, alpha, beta);
    double laPlaceCorr = corr->laPlaceWaveFunction(r,beta);

    for(int i=0;i<nParticles;i++) {
        gradSlater = slater->gradientWaveFunction(r,i,1,alpha,beta);
        gradCorr = corr->gradientWaveFunction(r,i,beta);
        gradProduct += dot(gradSlater, gradCorr);
    }

    double kineticEnergy = -0.5*(laPlaceSlater+laPlaceCorr+2*gradProduct);

    return kineticEnergy;
}


double Hamiltonian::potentialEnergy(const mat &r)
{
    // Potential energy (1/r part)
    double potentialEnergy = 0;
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j)*r(i,j);
        }
        potentialEnergy -= charge / sqrt(rSingleParticle);
    }

    // Contribution from electron-electron potential (1/rij part)
    double r12 = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < i; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += pow((r(i,k) - r(j,k)),2);
            }
            potentialEnergy += 1 / sqrt(r12);
        }
    }


    return potentialEnergy;
}

double Hamiltonian::analyticEnergyHe(const mat &r, const double &alpha, const double &beta)
{
    double r1 = sqrt(norm(r.row(0),2));
    double r2 = sqrt(norm(r.row(1),2));
    double r12 = norm((r.row(1) - r.row(0)),2);
    double dot_r12;
    for (int i = 0; i < nDimensions; i++) { dot_r12 += r(0,i)*r(1,i); }
    double energy_l1 = (alpha-charge)*(1/r1 + 1/r2) + 1/r12 - pow(alpha,2);
    double energy_part = ((alpha*(r1+r2))/r12)*(1-dot_r12/(r1*r2)) - 1/(2*pow((1+beta*r12),2)) - 2/r12 + (2*beta)/(1+beta*r12);
    double energy_l2 = energy_l1 + 1/(2*pow((1+beta*r12),2))*energy_part;

    return energy_l2;
}
