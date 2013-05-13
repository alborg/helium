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

//dPsi/dalpha or dPsi/dbeta
double Hamiltonian::dPsi(int var, const mat &r, const double &alpha, const double &beta, slaterDeterminant *slater, correlation *corr) {

    double dPsi = 0;

    if(var == 1) { //derivative with respect to alpha:
        dPsi = slater->dWaveFunction_dalpha(r,alpha);
    }

    return dPsi;
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

    double kineticEnergy = -0.5*(laPlaceSlater);//+laPlaceCorr+2*gradProduct);

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
//    double r12 = 0;
//    for(int i = 0; i < nParticles; i++) {
//        for(int j = 0; j < i; j++) {
//            r12 = 0;
//            for(int k = 0; k < nDimensions; k++) {
//                r12 += pow((r(i,k) - r(j,k)),2);
//            }
//            potentialEnergy += 1 / sqrt(r12);
//        }
//    }


    return potentialEnergy;
}

double Hamiltonian::analyticEnergyHe(const mat &r, const double &alpha, const double &beta)
{
    vec r1vec = zeros<vec>(nDimensions);
    vec r2vec = zeros<vec>(nDimensions);
    double r1 = 0;
    double r2 = 0;
    double r12 = 0;
     double r12dot = 0;
    double Energy = 0;

    for(int d=0; d<nDimensions; d++) {
        r1vec(d) = r(0,d);
        r2vec(d) = r(1,d);
        r1 += r(0,d)*r(0,d);
        r2 += r(1,d)*r(1,d);
        r12 += pow((r1vec(d) - r2vec(d)),2);
         r12dot += r(0,d)*r(1,d);
    }
    r12 = sqrt(r12);
    r1 = sqrt(r1);
    r2 = sqrt(r2);

    double EL1 = (alpha - charge)*(1/r1 + 1/r2) + 1/r12 - alpha*alpha;

    //With corr terms
    Energy = -pow(alpha, 2) + (alpha - charge)*(1.0/r2 + 1.0/r1) + (alpha*(1 - r12dot/(r1*r2))*(r1 + r2)/r12 + 2*beta/(beta*r12 + 1) - 1/((beta*r12 + 1)*(2*beta*r12 + 2)) - 2/r12)/((beta*r12 + 1)*(2*beta*r12 + 2)) + 1.0/r12;
    //Without corr
    Energy = alpha*alpha - 2*alpha*(charge - 5/16);

    return Energy;
}

double Hamiltonian::analyticdEnergyHe(const mat &r, const double &alpha, const double &beta)
{
    double r1 = 0;
    double r2 = 0;
    double r12 = 0;
    double r12dot = 0;

    for(int d=0; d<nDimensions; d++) {
        r1 += r(0,d)*r(0,d);
        r2 += r(1,d)*r(1,d);
        r12 += pow((r(0,d) - r(1,d)),2);
        r12dot += r(0,d)*r(1,d);
    }
    r12 = sqrt(r12);
    r1 = sqrt(r1);
    r2 = sqrt(r2);

    //With corr terms
    double dE = -2*alpha + 1.0/r2 + (1 - r12dot/(r1*r2))*(r1 + r2)/(r12*(beta*r12 + 1)*(2*beta*r12 + 2)) + 1.0/r1;

    //Without corr terms
    dE = 2*alpha - 2*(charge - 5/16);

    return dE;
}
