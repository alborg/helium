#include "WaveFunction.h"
#include "lib.h"

#include <armadillo>


using namespace arma;
using namespace std;

WaveFunction::WaveFunction(int &nParticles_, int &nDimensions_) :

    nDimensions(nDimensions_),
    nParticles(nParticles_)


{
}



double WaveFunction::gradientWaveFunction(const mat &r, int particle, int dimension, double alpha, double beta) {

    double variable = r(particle,dimension);
    double rtot = 0;
    double derivate = 0;

    for (int j=0; j<nDimensions; j++) { rtot += r(particle,j)*r(particle,j); }

    if(particle == 0 || particle == nParticles/2) {
        derivate = dPsi1s(rtot, variable, alpha); //n=1,l=0,ml=0
    }
    if(particle == 1 || particle == 1+nParticles/2) {
        derivate = dPsi2s(rtot, variable, alpha);//n=2,l=0,ml=0
    }
    if(particle == 2 || particle == 2+nParticles/2) {
        derivate = r(particle,1)*dPsi2p(rtot, variable, alpha);//n=2,l=1,ml=-1
        if(dimension == 1) derivate += psi2p(rtot, alpha);
    }
    if(particle == 3 || particle == 3+nParticles/2) {
        derivate = r(particle,0)*dPsi2p(rtot, variable, alpha);//n=2,l=1,ml=0
        if(dimension == 0) derivate += psi2p(rtot, alpha);
    }
    if(particle == 4 || particle == 4+nParticles/2) {
        derivate = r(particle,2)*dPsi2p(rtot, variable, alpha);//n=2,l=1,ml=1
        if(dimension == 2) derivate += psi2p(rtot, alpha);
    }

    return derivate;

}


//Wavefunction, 1s state
double WaveFunction::psi1s(double r, double alpha) {

    return exp(-alpha*r);

}

//First derivative of wavefunction, 1s state
double WaveFunction::dPsi1s(double rtot, double variable, double alpha) {

   return -(alpha*variable/rtot)*exp(-alpha*rtot);

}



//Wavefunction, 2s state
double WaveFunction::psi2s(double r, double alpha) {

    return (1-(alpha*r)/2)*exp(-alpha*r/2);

}

//First derivative of wavefunction, 2s state
double WaveFunction::dPsi2s(double rtot, double variable, double alpha) {

    return (alpha*variable/rtot)*(alpha*rtot/4 - 1)*exp(-alpha*rtot/2);

}

//Wavefunction, 2p state
double WaveFunction::psi2p(double r, double alpha) {

    return alpha*r*exp(-alpha*r/2);

}


//First derivative of wavefunction, 2p state
double WaveFunction::dPsi2p(double rtot, double variable, double alpha) {

    return (alpha*variable/rtot)*(1-rtot/2)*exp(-alpha*rtot/2);

}

double WaveFunction::jastrowFactor(const mat &r, double beta) {

    rowvec r12;
    double r12norm = 0;
    double jastrow = 0;

    for(int k=0;k<nParticles;k++) {
        for(int l=0;l<nParticles;l++) {
            if(k<l) {
                r12 = r.row(k) - r.row(l);
                r12norm = 0;
                for(int j = 0; j < nDimensions; j++) {
                    r12norm +=  r12(j)*r12(j);
                }
                jastrow += r12norm / (2 * (1 + beta * r12norm));
            }
        }
    }


    return exp(jastrow);

}
