#include "WaveFunction.h"
#include "lib.h"

#include <armadillo>


using namespace arma;
using namespace std;

WaveFunction::WaveFunction(int &nDimensions_, int &nProtons_, int &nElectrons_) :

    nDimensions(nDimensions_),
    nProtons(nProtons_),
    nElectrons(nElectrons_),
    nParticles(nProtons*nElectrons)

{
}


double WaveFunction::waveFunction(const mat &r, const mat &rProtons, double alpha, double beta) {

    double rp = 0;
    double expP = 0;
    double wave = 1;

    for(int e=0; e<nParticles; e++) {
        expP = 0;
        for(int p=0; p<nProtons; p++) {
            rp = 0;
            for(int d=0; d<nDimensions; d++) rp += (r(e,d) - rProtons(p,d))*(r(e,d) - rProtons(p,d));
            expP += exp(-alpha*sqrt(rp));
        }
        wave = wave * expP;
    }
    cout <<"wave: "<<wave<<endl;

    double jastrow = jastrowFactor(r,beta);

    cout<<"jastrow: "<<jastrow<<endl;

    return wave*jastrow;
}



double WaveFunction::jastrowFactor(const mat &r, double beta) {

    rowvec r12;
    double r12norm = 0;
    double jastrow = 0;
    double a = 1;
    vec spins = ones<vec>(nParticles);

    if(nElectrons > 1) {
        for(int e=0; e<nElectrons/2; e++) {
            spins(e) = 1;
            spins(nElectrons/2 + e) = -1;
            spins(nElectrons + e) = 1;
            spins(nElectrons + nElectrons/2 + e) = -1;
        }
    }

    for(int k=1;k<nParticles;k++) {
        for(int l=0;l<k;l++) {
            if(nElectrons > 1) {
                if(spins(k)*spins(l) > 0) a = 0.25;
                if(spins(k)*spins(l) < 0) a = 0.5;
            }
            r12 = r.row(k) - r.row(l);
            r12norm = 0;
            for(int d = 0; d < nDimensions; d++) r12norm +=  r12(d)*r12(d);
            jastrow += a * sqrt(r12norm) / (2 * (1 + beta * sqrt(r12norm)));

        }
    }



    return exp(jastrow);

}
