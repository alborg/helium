#include "WaveFunction.h"
#include "lib.h"
#include "slaterdeterminant.h"

#include <armadillo>


using namespace arma;
using namespace std;

WaveFunction::WaveFunction(int &nParticles_, int &nDimensions_) :

    nDimensions(nDimensions_),
    nParticles(nParticles_),
    slater(new slaterDeterminant(nParticles,nDimensions))


{
}

void WaveFunction::buildDeterminant(const mat &r, double alpha) {

    //Build the Slater determinant
    slater->buildDeterminant(r, alpha);

}

double WaveFunction::waveFunction(const mat &r, double alpha_, double beta_) {
    alpha = alpha_;
    beta = beta_;
    double argument = 0;
    double waveFunc = 0;

    //Helium
    if(nParticles==2) {
//        for(int i = 0; i < nParticles; i++) {
//            double rSingleParticle = 0;
//            for(int j = 0; j < nDimensions; j++) {
//                rSingleParticle += r(i,j) * r(i,j);
//            }
//            argument += sqrt(rSingleParticle);
//        }
//        waveFunc = exp(-argument * alpha) * jastrowFactor(r);

        waveFunc = slater->getDeterminant(r, alpha);//* jastrowFactor(r);
    }

//    Beryllium
   else if(nParticles ==4) {
        //waveFunc = slater->beryllium(r, alpha); //* jastrowFactor(r);
        //cout << waveFunc << endl;
        waveFunc = slater->getDeterminant(r, alpha);//* jastrowFactor(r);
        //cout << waveFunc << endl;
        //cout << "------" << endl;
    }

   else waveFunc = slater->getDeterminant(r, alpha); //* jastrowFactor(r);
   // cout << waveFunc << endl;

     return waveFunc;
}


//Wavefunction, 1s state
double WaveFunction::psi1s(double r, double alpha) {

    return exp(-alpha*r);

}


//Wavefunction, 2s state
double WaveFunction::psi2s(double r, double alpha) {

    return (1-(alpha*r)/2)*exp(-alpha*r/2);

}

//Wavefunction, 2p state
double WaveFunction::psi2p(double r, double alpha) {

    return alpha*r*exp(-alpha*r/2);

}


double WaveFunction::jastrowFactor(const mat &r) {

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
