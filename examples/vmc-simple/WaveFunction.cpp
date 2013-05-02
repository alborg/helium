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



//vec WaveFunction::gradientWaveFunction(const mat &r, int i, double alpha, double beta) {

//    double rtot = 0;
//    vec derivate = zeros<vec>(3);

//    for (int j=0; j<nDimensions; j++) { rtot += r(i,j)*r(i,j); }

//    if(i == 0 || i == nParticles/2) {
//        derivate = dPsi1s(rtot, i, r, alpha); //n=1,l=0,ml=0
//    }
//    if(i == 1 || i == 1+nParticles/2) {
//        derivate = dPsi2s(rtot, i, r, alpha);//n=2,l=0,ml=0
//    }
//    if(i == 2 || i == 2+nParticles/2) {
//        derivate = dPsi2p_1(rtot, i, r, alpha);//n=2,l=1,ml=-1
//    }
//    if(i == 3 || i == 3+nParticles/2) {
//        derivate = dPsi2p0(rtot, i, r, alpha);//n=2,l=1,ml=0
//    }
//    if(i == 4 || i == 4+nParticles/2) {
//        derivate = dPsi2p1(rtot, i, r, alpha);//n=2,l=1,ml=1
//    }

//    return derivate;

//}



//double WaveFunction::laPlaceWaveFunction(const mat &r, int i, double alpha, double beta) {

//    double rtot = 0;
//    double derivate = 0;

//    for (int j=0; j<nDimensions; j++) { rtot += r(i,j)*r(i,j); }

//    if(i == 0 || i == nParticles/2) {
//        derivate = d2Psi1s(rtot, i, r, alpha); //n=1,l=0,ml=0
//    }
//    if(i == 1 || i == 1+nParticles/2) {
//        derivate = d2Psi2s(rtot, i, r, alpha);//n=2,l=0,ml=0
//    }
//    if(i == 2 || i == 2+nParticles/2) {
//        derivate = d2Psi2p_1(rtot, i, r, alpha);//n=2,l=1,ml=-1
//    }
//    if(i == 3 || i == 3+nParticles/2) {
//        derivate = d2Psi2p0(rtot, i, r, alpha);//n=2,l=1,ml=0
//    }
//    if(i == 4 || i == 4+nParticles/2) {
//        derivate = d2Psi2p1(rtot, i, r, alpha);//n=2,l=1,ml=1
//    }

//    return derivate;

//}



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

    double der = (alpha/rtot)*(alpha*rtot - 2)*exp(-alpha*rtot);

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


//Wavefunction, 2p ml=0 state
double WaveFunction::psi2p0(double rtot, int i, const mat &r, double alpha) {

    double psi = r(i,0)*alpha*rtot*exp(-alpha*rtot/2);
    return psi;

}


//Wavefunction, 2p ml=-1 state
double WaveFunction::psi2p_1(double rtot, int i, const mat &r, double alpha) {

    double psi = r(i,1)*alpha*rtot*exp(-alpha*rtot/2);
    return psi;
}

//Wavefunction, 2p ml=1 state
double WaveFunction::psi2p1(double rtot, int i, const mat &r, double alpha) {

    double psi = r(i,2)*alpha*rtot*exp(-alpha*rtot/2);
    return psi;
}


//First derivative of wavefunction, 2p ml=0 state
vec WaveFunction::dPsi2p0(double rtot, int i, const mat &r, double alpha) {

    double x=r(i,0);
    double y=r(i,1);
    double z=r(i,2);

    vec der = zeros<vec>(3,1);

    der(0) = -alpha*(alpha*rtot*pow(x, 2) - 2*pow(rtot, 2) - 2*pow(x, 2))*exp(-alpha*rtot/2)/(2*rtot);
    der(1) = -alpha*x*y*(alpha*rtot - 2)*exp(-alpha*rtot/2)/(2*rtot);
    der(2) = -alpha*x*z*(alpha*rtot - 2)*exp(-alpha*rtot/2)/(2*rtot);


    return der;
}


//First derivative of wavefunction, 2p ml=-1 state
vec WaveFunction::dPsi2p_1(double rtot, int i, const mat &r, double alpha) {

    double x=r(i,0);
    double y=r(i,1);
    double z=r(i,2);

    vec der = zeros<vec>(3,1);

    der(0) = -alpha*x*y*(alpha*rtot - 2)*exp(-alpha*rtot/2)/(2*rtot);
    der(1) = -alpha*(alpha*rtot*pow(y, 2) - 2*pow(rtot, 2) - 2*pow(y, 2))*exp(-alpha*rtot/2)/(2*rtot);
    der(2) = -alpha*y*z*(alpha*rtot - 2)*exp(-alpha*rtot/2)/(2*rtot);

    return der;
}


//First derivative of wavefunction, 2p ml=1 state
vec WaveFunction::dPsi2p1(double rtot, int i, const mat &r, double alpha) {

    double x=r(i,0);
    double y=r(i,1);
    double z=r(i,2);

    vec der = zeros<vec>(3,1);

    der(0) = -alpha*x*z*(alpha*rtot - 2)*exp(-alpha*rtot/2)/(2*rtot);
    der(1) = -alpha*y*z*(alpha*rtot - 2)*exp(-alpha*rtot/2)/(2*rtot);
    der(2) = -alpha*(alpha*rtot*pow(z, 2) - 2*pow(rtot, 2) - 2*pow(z, 2))*exp(-alpha*rtot/2)/(2*rtot);

    return der;
}

//Second derivative of wavefunction, 2p ml=0 state
double WaveFunction::d2Psi2p0(double rtot, int i, const mat &r, double alpha) {
    double x=r(i,0);
    double y=r(i,1);
    double z=r(i,2);

    vec der2 = zeros<vec>(3,1);

    der2(0) = alpha*x*(pow(alpha, 2)*pow(rtot, 2)*pow(x, 2) - 6*alpha*pow(rtot, 3) - 2*alpha*rtot*pow(x, 2) + 12*pow(rtot, 2) - 4*pow(x, 2))*exp(-alpha*rtot/2)/(4*pow(rtot, 3));
    der2(1) = alpha*x*(pow(alpha, 2)*pow(rtot, 2)*pow(y, 2) - 2*alpha*pow(rtot, 3) - 2*alpha*rtot*pow(y, 2) + 4*pow(rtot, 2) - 4*pow(y, 2))*exp(-alpha*rtot/2)/(4*pow(rtot, 3));
    der2(2) = alpha*x*(pow(alpha, 2)*pow(rtot, 2)*pow(z, 2) - 2*alpha*pow(rtot, 3) - 2*alpha*rtot*pow(z, 2) + 4*pow(rtot, 2) - 4*pow(z, 2))*exp(-alpha*rtot/2)/(4*pow(rtot, 3));

   return sum(der2);

}


//Second derivative of wavefunction, 2p ml=-1 state
double WaveFunction::d2Psi2p_1(double rtot, int i, const mat &r, double alpha) {
    double x=r(i,0);
    double y=r(i,1);
    double z=r(i,2);
    vec der = zeros<vec>(3,1);
    der(0) = alpha*y*(pow(alpha, 2)*pow(rtot, 2)*pow(x, 2) - 2*alpha*pow(rtot, 3) - 2*alpha*rtot*pow(x, 2) + 4*pow(rtot, 2) - 4*pow(x, 2))*exp(-alpha*rtot/2)/(4*pow(rtot, 3));
    der(1) = alpha*y*(pow(alpha, 2)*pow(rtot, 2)*pow(y, 2) - 6*alpha*pow(rtot, 3) - 2*alpha*rtot*pow(y, 2) + 12*pow(rtot, 2) - 4*pow(y, 2))*exp(-alpha*rtot/2)/(4*pow(rtot, 3));
    der(2) = alpha*y*(pow(alpha, 2)*pow(rtot, 2)*pow(z, 2) - 2*alpha*pow(rtot, 3) - 2*alpha*rtot*pow(z, 2) + 4*pow(rtot, 2) - 4*pow(z, 2))*exp(-alpha*rtot/2)/(4*pow(rtot, 3));

   return sum(der);

}

//Second derivative of wavefunction, 2p ml=1 state
double WaveFunction::d2Psi2p1(double rtot, int i, const mat &r, double alpha) {
    double x=r(i,0);
    double y=r(i,1);
    double z=r(i,2);
    vec der = zeros<vec>(3,1);
    der(0) = alpha*z*(pow(alpha, 2)*pow(rtot, 2)*pow(x, 2) - 2*alpha*pow(rtot, 3) - 2*alpha*rtot*pow(x, 2) + 4*pow(rtot, 2) - 4*pow(x, 2))*exp(-alpha*rtot/2)/(4*pow(rtot, 3));
    der(1) = alpha*z*(pow(alpha, 2)*pow(rtot, 2)*pow(y, 2) - 2*alpha*pow(rtot, 3) - 2*alpha*rtot*pow(y, 2) + 4*pow(rtot, 2) - 4*pow(y, 2))*exp(-alpha*rtot/2)/(4*pow(rtot, 3));
    der(2) = alpha*z*(pow(alpha, 2)*pow(rtot, 2)*pow(z, 2) - 6*alpha*pow(rtot, 3) - 2*alpha*rtot*pow(z, 2) + 12*pow(rtot, 2) - 4*pow(z, 2))*exp(-alpha*rtot/2)/(4*pow(rtot, 3));

   return sum(der);

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
