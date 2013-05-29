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



//Wavefunction, 1s state
double WaveFunction::psi1s(double rtot, double alpha) {

    double psi = exp(-alpha*rtot);
    return psi;

}

//First derivative of wavefunction, 1s state
vec WaveFunction::dPsi1s(double rtot, int i, const mat &r, double alpha) {

    vec der = zeros<vec>(3,1);
    der(0) = -alpha*r(i,0)*exp(-alpha*rtot)/rtot;
    der(1) = -alpha*r(i,1)*exp(-alpha*rtot)/rtot;
    der(2) = -alpha*r(i,2)*exp(-alpha*rtot)/rtot;

   return der;

}



//Second derivative of wavefunction, 1s state
double WaveFunction::d2Psi1s(double rtot, double alpha) {

    double der = alpha*(alpha*rtot - 2)*exp(-alpha*rtot)/rtot;

   return der;

}


//Wavefunction, 2s state
double WaveFunction::psi2s(double rtot, double alpha) {

    double psi = (1-(alpha*rtot)/2)*exp(-alpha*rtot/2);
    return psi;

}

//First derivative of wavefunction, 2s state
vec WaveFunction::dPsi2s(double rtot, int i, const mat &r, double alpha) {

    vec der = zeros<vec>(3,1);

    der(0) = alpha*r(i,0)*(alpha*rtot - 4)*exp(-alpha*rtot/2)/(4*rtot);
    der(1) = alpha*r(i,1)*(alpha*rtot - 4)*exp(-alpha*rtot/2)/(4*rtot);
    der(2) = alpha*r(i,2)*(alpha*rtot - 4)*exp(-alpha*rtot/2)/(4*rtot);

    return der;
}



//Second derivative of wavefunction, 2s state
double WaveFunction::d2Psi2s(double rtot, double alpha) {

    double der = -(alpha/(8*rtot))*(alpha*rtot-8)*(alpha*rtot-2)*exp(-alpha*rtot/2);

    return der;

}

