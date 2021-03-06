#include "WaveFunction.h"
#include "lib.h"
#include <armadillo>


using namespace arma;
using namespace std;

WaveFunction::WaveFunction(int &nParticles_, int &nDimensions_) :

    nDimensions(nDimensions_),
    nParticles(nParticles_),
    slater(new slaterDeterminant(nParticles,nDimensions))


{
}

//Compute wavefunction of He:
double WaveFunction::waveFunction(const mat &r, double alpha_, double beta_) {
    alpha = alpha_;
    beta = beta_;
    double argument = 0;
    double waveFunc = 0;

    for(int i = 0; i < nParticles; i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        argument += sqrt(rSingleParticle);
    }
    waveFunc = exp(-argument * alpha) * jastrowFactor(r); //Both parts of wavefunction




    return waveFunc;
}

//Correlation part of wavefunction:
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
                r12norm = sqrt(r12norm);
                jastrow += r12norm / (2 * (1 + beta * r12norm));
            }
        }
    }


    return exp(jastrow);

}
