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


//Find the local energy (expectation value of the energy)
double Hamiltonian::localEnergy(const mat &r, const double &alpha, const double &beta, slaterDeterminant *slater, correlation *corr)
{


    double kinEnergy = kineticEnergy(r, alpha, beta, slater,corr);
    double potEnergy = potentialEnergy(r);

    return kinEnergy + potEnergy;
}

//dPsi/dalpha/Psi and dPsi/dbeta/Psi
vec Hamiltonian::dPsi(const mat &r, double alpha, double beta, slaterDeterminant *slater, correlation *corr) {

    vec dPsi = zeros<vec>(2,1);

    //First derivative of wavefunction wrt alpha
    slater->buildDeterminant(r,alpha,beta);
    double wf = slater->getDeterminant();

    double alphaPlus, alphaMinus;
    alphaPlus = alphaMinus = alpha;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    alphaPlus = alpha+h;
    alphaMinus = alpha-h;
    slater->buildDeterminant(r,alphaMinus,beta);
    waveFunctionMinus = slater->getDeterminant();
    slater->buildDeterminant(r,alphaPlus,beta);
    waveFunctionPlus = slater->getDeterminant();
    dPsi(0) = (waveFunctionPlus - waveFunctionMinus)/(2*wf*h);

    slater->buildDeterminant(r,alpha,beta);

    //First derivative of wavefunction wrt beta
    wf = corr->jastrow(r,beta);
    double betaPlus, betaMinus;
    betaPlus = betaMinus = beta;

    betaPlus = beta+h;
    betaMinus = beta-h;
    waveFunctionMinus = corr->jastrow(r,betaMinus);
    waveFunctionPlus = corr->jastrow(r,betaPlus);
    dPsi(1) = (waveFunctionPlus - waveFunctionMinus)/(2*wf*h);



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

