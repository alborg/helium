#include "hamiltonian.h"
#include "slaterdeterminant.h"
#include "correlation.h"

Hamiltonian::Hamiltonian(int nDimensions_, double h_, double h2_, int charge_, int nProtons_, int nElectrons_, double R_, const mat &rProtons_) :

    nDimensions(nDimensions_),
    nProtons(nProtons_),
    nElectrons(nElectrons_),
    nParticles(nProtons*nElectrons),
    h(h_),
    h2(h2_),
    charge(charge_),
    R(R_),
    rProtons(rProtons_)
{
}


//Find the local energy (expectation value of the energy) numerically
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
    //slater->buildDeterminant(r,alpha);
    double wf = slater->getDeterminant();

    double alphaPlus, alphaMinus;
    alphaPlus = alphaMinus = alpha;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    alphaPlus = alpha+h;
    alphaMinus = alpha-h;
    slater->buildDeterminant(r,alphaMinus);
    waveFunctionMinus = slater->getDeterminant();
    slater->buildDeterminant(r,alphaPlus);
    waveFunctionPlus = slater->getDeterminant();
    dPsi(0) = (waveFunctionPlus - waveFunctionMinus)/(2*wf*h);

    slater->buildDeterminant(r,alpha);

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

    vec gradSlater1 = zeros<vec>(nDimensions,1);
    vec gradSlater2 = zeros<vec>(nDimensions,1);
    vec gradCorr = zeros<vec>(nDimensions,1);
    double laPlaceSlater1 = 0;
    double laPlaceSlater2 = 0;
    double gradProduct = 0;


    double laPlaceCorr = corr->laPlaceWaveFunction(r,beta);

    for(int i=0;i<nParticles;i++) gradCorr = corr->gradientWaveFunction(r,i,beta);

    for(int i=0;i<nElectrons;i++) { //Atom1
        laPlaceSlater1 += slater->laPlaceWaveFunction(r, i, alpha);
        gradSlater1 += slater->gradientWaveFunction(r,i,1,alpha);
    }

    for(int j=nElectrons;j<nParticles;j++) { //Atom2
        laPlaceSlater2 += slater->laPlaceWaveFunction(r, j, alpha);
        gradSlater2 += slater->gradientWaveFunction(r,j,1,alpha);
    }

    gradProduct = dot(gradSlater1, gradCorr) + dot(gradSlater2, gradCorr) + dot(gradSlater1, gradSlater2);

    double kineticEnergy = -0.5*(laPlaceSlater1+laPlaceSlater2+laPlaceCorr+2*gradProduct);

    return kineticEnergy;
}


double Hamiltonian::potentialEnergy(const mat &r)
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
