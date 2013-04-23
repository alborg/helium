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
        for(int g = 0; g < nDimensions; g++) {
            rSingleParticle += r(i,g) * r(i,g);
        }
        rs(i) = sqrt(rSingleParticle);
    }

    //Make Slater determinants (spin up and spin down)
    for(int j=0; j<nParticles/2; j++) { //Rows: State
        for(int i=0; i<nParticles/2; i++) { //Cols: Particle (position)
            if (j == 0) { //n=1,l=0,ml=0
                slaterMatrixUp(j,i) = function->psi1s(rs[i], alpha);
                slaterMatrixDown(j,i) = function->psi1s(rs[nParticles/2+i], alpha);
            }
            if (j == 1) {//n=2,l=0,ml=0
                slaterMatrixUp(j,i) = function->psi2s(rs[i], alpha);
                slaterMatrixDown(j,i) = function->psi2s(rs[nParticles/2+i], alpha);
            }
            if (j == 2) { //n=2,l=1,ml=-1
                slaterMatrixUp(j,i) = function->psi2p_1(rs[i], i, r, alpha);
                slaterMatrixDown(j,i) = function->psi2p_1(rs[nParticles/2+i], nParticles/2+i, r, alpha);
            }
            if (j == 3) { //n=2,l=1,ml=0
                slaterMatrixUp(j,i) = function->psi2p0(rs[i], i, r, alpha);
                slaterMatrixDown(j,i) = function->psi2p0(rs[nParticles/2+i], nParticles/2+i, r, alpha);
            }
            if (j == 4) { //n=2,l=1,ml=1
                slaterMatrixUp(j,i) = function->psi2p1(rs[i], i, r, alpha);
                slaterMatrixDown(j,i) = function->psi2p1(rs[nParticles/2+i], nParticles/2+i, r, alpha);
            }
        }
    }

    invSlaterMatrixUp = inv(slaterMatrixUp);
    invSlaterMatrixDown = inv(slaterMatrixDown);


}


double slaterDeterminant::getDeterminant() {

    return det(slaterMatrixUp)*det(slaterMatrixDown);

}



double slaterDeterminant::getRatioDeterminant(int i, const mat &r, double alpha, double beta) {

    vec ratio = zeros<vec>(nParticles/2,1);
    double rSingleParticle = 0;
    //Get rtot for particle's new position:
    for(int g = 0; g < nDimensions; g++) {
        rSingleParticle += r(i,g) * r(i,g);
    }
    double rtot = sqrt(rSingleParticle);

    vec updatedStates = getStates(r, i, rtot, alpha, beta); //Get states with new r for particle i

    if(i<nParticles/2) {//Particle spin up
        for(int j=0; j<nParticles/2; j++) { //States
            ratio(j) = updatedStates(j) * invSlaterMatrixUp(j,i);
        }
    }
    else { //Particle spin down
        for(int j=0; j<nParticles/2; j++) { //States
            ratio(j) = updatedStates(j) * invSlaterMatrixDown(j,i-nParticles/2);
        }
    }


    return sum(ratio);
}


vec slaterDeterminant::getStates(const mat &r, int i, double rtot, double alpha, double beta) {

    vec updatedStates = zeros<vec>(nParticles/2,1);

    updatedStates(0) = function->psi1s(rtot, alpha); //n=1,l=0,ml=0
    if(nParticles/2>1) updatedStates(1) = function->psi2s(rtot, alpha); //n=2,l=0,ml=0
    if(nParticles/2>2) updatedStates(2) = function->psi2p_1(rtot, i, r, alpha); //n=2,l=1,ml=-1
    if(nParticles/2>3) updatedStates(3) = function->psi2p0(rtot, i, r, alpha); //n=2,l=1,ml=0
    if(nParticles/2>4) updatedStates(4) = function->psi2p1(rtot, i, r, alpha); //n=2,l=1,ml=1

    return updatedStates;
}



void slaterDeterminant::updateDeterminant(const mat &rNew, const mat &rOld, int i, double &alpha_, double &beta_, double ratio) {


    double rtot = 0;
    //Get rtot for particles' position:
    for(int d = 0; d < nDimensions; d++) {
        rtot += rNew(i,d) * rNew(i,d);
    }
    rtot = sqrt(rtot);

    vec newStates = getStates(rNew, i, rtot, alpha_, beta_);
    vec sumSjUp = zeros<vec>(nParticles/2); //Sum over states(l)*d_lj for particles j
    vec sumSjDown = zeros<vec>(nParticles/2); //Sum over states(l)*d_lj for particles j

    for(int j=0; j<nParticles/2; j++) { //particle
        for(int l=0; l<nParticles/2; l++) { //state
            sumSjUp(j) += newStates(l) * invSlaterMatrixUp(l,j);
            sumSjDown(j) += newStates(l) * invSlaterMatrixDown(l,j);
        }
    }

    int particle = i;
    if(i>=nParticles/2) particle = i-nParticles/2;

    //Update inverse matrices:

    //All columns except column corresponding to particle i:
    if(i<nParticles/2) { //If particle i has spin up
        for(int k=0; k<nParticles/2; k++) { //Rows, inv matrix
            for (int j=0; j<nParticles/2; j++) { //Cols, inv matrix
                if(j != particle) invSlaterMatrixUp(k,j) = invSlaterMatrixUp(k,j) - (sumSjUp(j)/ratio)*invSlaterMatrixUp(k,particle);
                invSlaterMatrixDown(k,j) = invSlaterMatrixDown(k,j) - (sumSjDown(j)/ratio)*invSlaterMatrixDown(k,particle);
            }

        }
    }
    else {  //If particle i has spin down
        for(int k=0; k<nParticles/2; k++) { //Rows, inv matrix
            for (int j=0; j<nParticles/2; j++) { //Cols, inv matrix
                invSlaterMatrixUp(k,j) = invSlaterMatrixUp(k,j) - (sumSjUp(j)/ratio)*invSlaterMatrixUp(k,particle);
                if(j != particle) invSlaterMatrixDown(k,j) = invSlaterMatrixDown(k,j) - (sumSjDown(j)/ratio)*invSlaterMatrixDown(k,particle);
            }

        }
    }


    //Update column corresponding to particle i:
    for(int k=0; k<nParticles/2; k++) { //States (rows)
        if(i<nParticles/2) { invSlaterMatrixUp(k,particle) = (1/ratio)*invSlaterMatrixUp(k,particle); }
        else {invSlaterMatrixDown(k,particle) = (1/ratio)*invSlaterMatrixDown(k,particle); }
    }


}

vec slaterDeterminant::gradientWaveFunction(const mat &r, int i, double ratio, double alpha, double beta) {

    double rtot = 0;
    vec derivate = zeros<vec>(nDimensions,1);
    vec invMatrix = zeros<vec>(nParticles/2,1);
    for(int p=0; p<nParticles/2; p++) {
        if(i<nParticles/2) { invMatrix(p) = invSlaterMatrixUp(p,i); }
        else { invMatrix(p) = invSlaterMatrixDown(p,i); }
    }

    for (int d=0; d<nDimensions; d++) { rtot += r(i,d)*r(i,d); }

    derivate += function->dPsi1s(rtot, i, r, alpha)*(1/ratio)*invMatrix(0); //n=1,l=0,ml=0
    if(nParticles > 1) {
        derivate += function->dPsi2s(rtot, i, r, alpha)*(1/ratio)*invMatrix(1);
    } //n=2,l=0,ml=0
    if(nParticles > 2) {
        derivate += function->dPsi2p_1(rtot, i, r, alpha)*(1/ratio)*invMatrix(2);
    } //n=2,l=1,ml=-1
    if(nParticles > 3) {
        derivate += function->dPsi2p0(rtot, i, r, alpha)*(1/ratio)*invMatrix(3);
    } //n=2,l=1,ml=0
    if(nParticles > 4) {
        derivate += function->dPsi2p1(rtot, i, r, alpha)*(1/ratio)*invMatrix(4);
    } //n=2,l=1,ml=1

    return derivate;

}


double slaterDeterminant::laPlaceWaveFunction(const mat &r, double alpha, double beta) {

    double kineticEnergy = 0;
    double rtot = 0;
    double derivate = 0;

    for(int i = 0; i < nParticles; i++) { //Particles

        for (int d=0; d<nDimensions; d++) { rtot += r(i,d)*r(i,d); } //Get r for particle

        for (int j=0; j<nParticles/2; j++) { //States

            //Get laplacian for each state:
            if(j == 0) {
                derivate = function->d2Psi1s(rtot, i, r, alpha); //n=1,l=0,ml=0
            }
            if(j == 1) {
                derivate = function->d2Psi2s(rtot, i, r, alpha);//n=2,l=0,ml=0
            }
            if(j == 2) {
                derivate = function->d2Psi2p_1(rtot, i, r, alpha);//n=2,l=1,ml=-1
            }
            if(j == 3) {
                derivate = function->d2Psi2p0(rtot, i, r, alpha);//n=2,l=1,ml=0
            }
            if(j == 4) {
                derivate = function->d2Psi2p1(rtot, i, r, alpha);//n=2,l=1,ml=1
            }

            if(i<nParticles/2) { derivate = derivate*invSlaterMatrixUp(j,i); }
            else { derivate = derivate*invSlaterMatrixDown(j,i-nParticles/2); }
            kineticEnergy += derivate;

        }
    }

    return kineticEnergy;

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


