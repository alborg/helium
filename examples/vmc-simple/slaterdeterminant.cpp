#include "slaterdeterminant.h"
#include "WaveFunction.h"

#include <armadillo>


using namespace arma;
using namespace std;

slaterDeterminant::slaterDeterminant(int nParticles_, int nDimensions_):

    nDimensions(nDimensions_),
    nParticles(nParticles_),
    slaterMatrixUp(zeros(nParticles/2,nParticles/2)),
    slaterMatrixDown(zeros(nParticles/2,nParticles/2)),
    slaterMatrixUpNew(zeros(nParticles/2,nParticles/2)),
    slaterMatrixDownNew(zeros(nParticles/2,nParticles/2)),
    invSlaterMatrixUp(zeros(nParticles/2,nParticles/2)),
    invSlaterMatrixDown(zeros(nParticles/2,nParticles/2)),
    function(new WaveFunction(nParticles,nDimensions))

{
}


void slaterDeterminant::buildDeterminant(const mat &r, double &alpha_, double &beta_) {

    vec rs = zeros<vec>(nParticles,1);
    alpha = alpha_; //Set class variable alpha

    //Find |r| for each electron:
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        rs(i) = sqrt(rSingleParticle);
    }


    //Make Slater determinants (spin up and spin down)
    for(int i=0; i<nParticles/2; i++) { //Rows: State
        for(int k=0; k<nParticles/2; k++) { //Cols: Particle (position)
            if (i == 0) { //n=1,l=0,ml=0
                slaterMatrixUp(i,k) = function->psi1s(rs[k], alpha);
                slaterMatrixDown(i,k) = function->psi1s(rs[nParticles/2+k], alpha);
            }
            if (i == 1) {//n=2,l=0,ml=0
                slaterMatrixUp(i,k) = function->psi2s(rs[k], alpha);
                slaterMatrixDown(i,k) = function->psi2s(rs[nParticles/2+k], alpha);
            }
            if (i == 2) { //n=2,l=1,ml=-1
                slaterMatrixUp(i,k) = r(k,1)*function->psi2p(rs[k], alpha);
                slaterMatrixDown(i,k) = r(k,1)*function->psi2p(rs[nParticles/2+k], alpha);
            }
            if (i == 3) { //n=2,l=1,ml=0
                slaterMatrixUp(i,k) = r(k,0)*function->psi2p(rs[k], alpha);
                slaterMatrixDown(i,k) = r(k,0)*function->psi2p(rs[nParticles/2+k], alpha);
            }
            if (i == 4) { //n=2,l=1,ml=1
                slaterMatrixUp(i,k) = r(k,2)*function->psi2p(rs[k], alpha);
                slaterMatrixDown(i,k) = r(k,2)*function->psi2p(rs[nParticles/2+k], alpha);
            }
            invSlaterMatrixUp(i,k) = 1/slaterMatrixUp(i,k);
            invSlaterMatrixDown(i,k) = 1/slaterMatrixDown(i,k);
        }


    }

}


double slaterDeterminant::getDeterminant() {

    return det(slaterMatrixUp)*det(slaterMatrixDown);

}

void slaterDeterminant::updateDeterminant(const mat &r, int i, double &alpha_, double &beta_) {

    for (int g=0; g<nParticles; g++) {
        if(g != i) {
            for(int j=0; j<nDimensions; j++) {



            }
        }
    }


}


double slaterDeterminant::getRatioDeterminant(int i, const double &r) {

    vec updatedStates = zeros<vec>(nParticles/2,1);

    double rSingleParticle;
    //Get rtot for particle's new position:
    for(int j = 0; j < nDimensions; j++) {
        rSingleParticle += r(i,j) * r(i,j);
    }
    double rtot = sqrt(rSingleParticle);

    updatedStates(0) = function->psi1s(rtot, alpha); //n=1,l=0,ml=0
    if(i>0) updatedStates(1) = function->psi2s(rtot, alpha); //n=2,l=0,ml=0
    if(i>1) updatedStates(2) = r(i,1)*function->psi2p(rtot, alpha); //n=2,l=1,ml=-1
    if(i>2) updatedStates(3) = r(i,0)*function->psi2p(rtot, alpha); //n=2,l=1,ml=0
    if(i>3) updatedStates(4) = r(i,2)*function->psi2p(rtot, alpha); //n=2,l=1,ml=1

    for(int g=0; g<nParticles/2; g++) {
        if(i<nParticles/2) { updatedStates(g) *= invSlaterMatrixUp(i,g); }
        else { updatedStates(g) *= invSlaterMatrixDown(i,g); }
    }

    return sum(updatedStates);
}



double slaterDeterminant::beryllium(const mat &r, double &alpha_)  {

    double rs[nParticles];
    alpha = alpha_; //Set class variable alpha

    //Find |r| for each electron:
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        rs[i] = sqrt(rSingleParticle);
    }

    //Slater determinant, Be
    double waveFunction = (function->psi1s(rs[0], alpha)*function->psi2s(rs[1], alpha) - function->psi1s(rs[1], alpha)*function->psi2s(rs[0], alpha))*
                          (function->psi1s(rs[2], alpha)*function->psi2s(rs[3], alpha) - function->psi1s(rs[3], alpha)*function->psi2s(rs[2], alpha));

    return waveFunction;

}


