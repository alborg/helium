#include "slaterdeterminant.h"

#include <armadillo>


using namespace arma;
using namespace std;

slaterDeterminant::slaterDeterminant(int nParticles_, int nDimensions_):

    nDimensions(nDimensions_),
    nParticles(nParticles_)

{
}


void slaterDeterminant::buildDeterminant(const mat &r, double &alpha_) {

    double rs[nParticles];
    alpha = alpha_; //Set class variable alpha

    slaterMatrixUp = zeros(nParticles/2,nParticles/2);
    slaterMatrixDown = zeros(nParticles/2,nParticles/2);

    //Find |r| for each electron:
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        rs[i] = sqrt(rSingleParticle);
    }

    //Make Slater determinants (spin up and spin down)
    for(int i=0; i<nParticles/2; i++) { //Rows: State
        for(int k=0; k<nParticles/2; k++) { //Cols: Particle (position)
            if (i == 0) { //n=1,l=0,ml=0
                slaterMatrixUp(i,k) = psi1s(rs[k]);
                slaterMatrixDown(i,nParticles/2+k) = psi1s(rs[nParticles/2+k]);
            }
            if (i == 1) {//n=2,l=0,ml=0
                slaterMatrixUp(i,k) = psi2s(rs[k]);
                slaterMatrixDown(i,nParticles/2+k) = psi2s(rs[nParticles/2+k]);
            }
            if (i == 2) { //n=2,l=1,ml=-1
                slaterMatrixUp(i,k) = r(k,1)*psi2p(rs[k]);
                slaterMatrixDown(i,nParticles/2+k) = r(k,1)*psi2p(rs[nParticles/2+k]);
            }
            if (i == 3) { //n=2,l=1,ml=0
                slaterMatrixUp(i,k) = r(k,0)*psi2p(rs[k]);
                slaterMatrixDown(i,nParticles/2+k) = r(k,0)*psi2p(rs[nParticles/2+k]);
            }
            if (i == 4) { //n=2,l=1,ml=1
                slaterMatrixUp(i,k) = r(k,2)*psi2p(rs[k]);
                slaterMatrixDown(i,nParticles/2+k) = r(k,2)*psi2p(rs[nParticles/2+k]);
            }
        }

    }

}


double slaterDeterminant::getDeterminant(const mat &r, double &alpha_) {

    alpha = alpha_; //Set class variable alpha

    return det(slaterMatrixUp*slaterMatrixDown);

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
    double waveFunction = (psi1s(rs[0])*psi2s(rs[1]) - psi1s(rs[1])*psi2s(rs[0]))*
                          (psi1s(rs[2])*psi2s(rs[3]) - psi1s(rs[3])*psi2s(rs[2]));

    return waveFunction;

}


//Wavefunction, 1s state
double slaterDeterminant::psi1s(double r) {

    return exp(-alpha*r);

}


//Wavefunction, 2s state
double slaterDeterminant::psi2s(double r) {

    return (1-(alpha*r)/2)*exp(-alpha*r/2);

}

//Wavefunction, 2p state
double slaterDeterminant::psi2p(double r) {

    return alpha*r*exp(-alpha*r/2);

}
