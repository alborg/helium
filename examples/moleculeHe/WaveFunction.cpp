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

//Get wavefunction
double WaveFunction::waveFunction(const mat &r, const mat &rProtons, double alpha, double beta) {
    double rp;
    double rp1 = 0;
    double rp2 = 0;
    double rp3 = 0;
    double rp4 = 0;
    double expP = 0;
    double wave = 1;

    for(int d=0; d<nDimensions; d++) {
        rp1 += pow(r(0,d) - rProtons(0,d),2);
        rp2 += pow(r(0,d) - rProtons(1,d),2);
        rp3 += pow(r(1,d) - rProtons(0,d),2);
        rp4 += pow(r(1,d) - rProtons(1,d),2);
    }
    wave = (exp(-alpha*sqrt(rp1)) + exp(-alpha*sqrt(rp2))) * (exp(-alpha*sqrt(rp3)) + exp(-alpha*sqrt(rp4)));
//wave = 1;
//    for(int e=0; e<nParticles; e++) {
//        expP = 0;
//        rp = 0;
//        for(int p=0; p<nProtons; p++) {
//            rp = 0;
//            for(int d=0; d<nDimensions; d++) rp += (r(e,d) - rProtons(p,d))*(r(e,d) - rProtons(p,d));
//            expP += exp(-alpha*sqrt(rp));
//        }
//        wave = wave * expP;
//    }
    //cout <<"wave: "<<wave<<endl;

    double jastrow = jastrowFactor(r,beta);

    //cout<<"jastrow: "<<jastrow<<endl;

    return wave*jastrow;
}


//Get correlation factor (Jastrow factor)
double WaveFunction::jastrowFactor(const mat &r, double beta) {

    rowvec r12;
    double r12norm = 0;
    double jastrow = 0;

    for(int k=1;k<nParticles;k++) {
        for(int l=0;l<k;l++) {
            r12 = r.row(k) - r.row(l);
            r12norm = 0;
            for(int d = 0; d < nDimensions; d++) r12norm +=  r12(d)*r12(d);
            jastrow += sqrt(r12norm) / (2 * (1 + beta * sqrt(r12norm)));

        }
    }



    return exp(jastrow);

}
