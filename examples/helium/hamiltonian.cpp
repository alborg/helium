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


//Find the local energy (expectation value of the energy) numerically
double Hamiltonian::localEnergy(const mat &r, const double &alpha, const double &beta, WaveFunction *function)
{
    double kinEnergy = kineticEnergy(r, alpha, beta, function);
    double potEnergy = potentialEnergy(r);

    return kinEnergy + potEnergy;
}

//Find the kinetic energy part of the local energy
double Hamiltonian::kineticEnergy(const mat &r, const double &alpha, const double &beta, WaveFunction *function)
{

    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    double waveFunctionCurrent = function->waveFunction(r, alpha, beta); //Find wavefunction for r

    //Second derivative (del^2):

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
        for(int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            potentialEnergy += 1 / sqrt(r12);
        }
    }

    return potentialEnergy;
}

//Find derivative of wavefunction wrt alpha, beta
vec Hamiltonian::dPsi(const mat &r, double alpha, double beta, WaveFunction *function) {

    vec dPsi = zeros<vec>(2,1);

    double wf = function->waveFunction(r, alpha, beta); //Find wavefunction for r

    //First derivative of wavefunction wrt alpha


    double alphaPlus, alphaMinus;
    alphaPlus = alphaMinus = alpha;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    alphaPlus = alpha+h;
    alphaMinus = alpha-h;

    waveFunctionMinus = function->waveFunction(r, alphaMinus, beta);
    waveFunctionPlus = function->waveFunction(r, alphaPlus, beta);
    dPsi(0) = (waveFunctionPlus - waveFunctionMinus)/(2*wf*h);


    //First derivative of wavefunction wrt beta
    double betaPlus, betaMinus;
    betaPlus = betaMinus = beta;

    betaPlus = beta+h;
    betaMinus = beta-h;
    waveFunctionMinus = function->waveFunction(r,alpha, betaMinus);
    waveFunctionPlus = function->waveFunction(r,alpha, betaPlus);
    dPsi(1) = (waveFunctionPlus - waveFunctionMinus)/(2*wf*h);



    return dPsi;
}




//The analytic energy of He
double Hamiltonian::analyticEnergyHe(const mat &r, const double &alpha, const double &beta)
{
    double r1 = 0;
    double r2 = 0;
    double r12 = 0;
    double dot_r12 = 0;
    for(int d=0; d<nDimensions; d++) {
        r1 += pow(r(0,d),2);
        r2 += pow(r(1,d),2);
        r12 += pow(r(0,d) - r(1,d),2);
        dot_r12 += r(0,d)*r(1,d);
    }
    r1 = sqrt(r1);
    r2 = sqrt(r2);
    r12 = sqrt(r12);

    double energy_l1 = (alpha-charge)*(1/r1 + 1/r2) + 1/r12 - pow(alpha,2);
    double energy_part = ((alpha*(r1+r2))/r12)*(1-dot_r12/(r1*r2)) - 1/(2*pow((1+beta*r12),2)) - 2/r12 + (2*beta)/(1+beta*r12);
    double energy_l2 = energy_l1 + 1/(2*pow((1+beta*r12),2))*energy_part;

    return energy_l2;
}
