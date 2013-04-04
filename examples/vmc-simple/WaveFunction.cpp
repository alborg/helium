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
        for(int i = 0; i < nParticles; i++) {
            double rSingleParticle = 0;
            for(int j = 0; j < nDimensions; j++) {
                rSingleParticle += r(i,j) * r(i,j);
            }
            argument += sqrt(rSingleParticle);
        }
        waveFunc = exp(-argument * alpha) * jastrowFactor(r);
    }

//    //Beryllium
   else if(nParticles ==4) {
        waveFunc = slater->beryllium(r, alpha); //* jastrowFactor(r);
        //waveFunc = slater->getDeterminant(r, alpha);//* jastrowFactor(r);

    }

   else slater->buildDeterminant(r, alpha); //* jastrowFactor(r);


     return waveFunc;
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
