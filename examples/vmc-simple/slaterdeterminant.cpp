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
    double alpha = alpha_;

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
    for(int i=0; i<nParticles/2; i++) { //Rows: Particle (position)
        for(int j=0; j<nParticles/2; j++) { //Cols: State
            if (j == 0) { //n=1,l=0,ml=0
                slaterMatrixUp(i,j) = function->psi1s(rs[i], alpha);
                slaterMatrixDown(i,j) = function->psi1s(rs[nParticles/2+i], alpha);
            }
            if (j == 1) {//n=2,l=0,ml=0
                slaterMatrixUp(i,j) = function->psi2s(rs[i], alpha);
                slaterMatrixDown(i,j) = function->psi2s(rs[nParticles/2+i], alpha);
            }
            if (j == 2) { //n=2,l=1,ml=-1
                slaterMatrixUp(i,j) = function->psi2p_1(rs[i], i, r, alpha);
                slaterMatrixDown(i,j) = function->psi2p_1(rs[nParticles/2+i], nParticles/2+i, r, alpha);
            }
            if (j == 3) { //n=2,l=1,ml=0
                slaterMatrixUp(i,j) = function->psi2p0(rs[i], i, r, alpha);
                slaterMatrixDown(i,j) = function->psi2p0(rs[nParticles/2+i], nParticles/2+i, r, alpha);
            }
            if (j == 4) { //n=2,l=1,ml=1
                slaterMatrixUp(i,j) = function->psi2p1(rs[i], i, r, alpha);
                slaterMatrixDown(i,j) = function->psi2p1(rs[nParticles/2+i], nParticles/2+i, r, alpha);
            }
        }
    }

    invSlaterMatrixUp = inv(slaterMatrixUp);
    invSlaterMatrixDown = inv(slaterMatrixDown);


}


double slaterDeterminant::getDeterminant() {

//    cout << "matrise"<<endl;
//    cout << slaterMatrixUp << endl;
//    cout << slaterMatrixDown << endl;
//    cout <<"----------"<<endl;


    return det(slaterMatrixUp)*det(slaterMatrixDown);

}


double slaterDeterminant::getInvDeterminant() {

//    cout << "inv matr"<<endl;
//    cout << inv(invSlaterMatrixUp) << endl;
//    cout << inv(invSlaterMatrixDown) << endl;
//    cout <<"----------"<<endl;

    return 1/(det(invSlaterMatrixUp)*det(invSlaterMatrixDown));

}


double slaterDeterminant::getRatioDeterminant(int i, const mat &r, double alpha, double beta) {

    double ratio = 0;
    double rSingleParticle = 0;
    //Get rtot for particle i's new position:
    for(int d = 0; d < nDimensions; d++) {
        rSingleParticle += r(i,d) * r(i,d);
    }
    double rtot = sqrt(rSingleParticle);

    vec updatedStates = getStates(r, i, rtot, alpha, beta); //Get states with new r for particle i

    if(i<nParticles/2) { //Particle spin up
        for(int j=0; j<nParticles/2; j++) { //States
            ratio += updatedStates(j) * invSlaterMatrixUp(j,i);
        }
    }
    else { //Particle spin down
        for(int j=0; j<nParticles/2; j++) { //States
            ratio += updatedStates(j) * invSlaterMatrixDown(j,i-nParticles/2);
        }
    }


    return ratio;
}


double slaterDeterminant::getRatioDeterminantNum(int i, const mat &rOld, const mat &rNew, double alpha, double beta) {

    return beryllium(rNew, alpha) / beryllium(rOld, alpha);

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

    int particle = i;
    if(i>=nParticles/2) particle = i-nParticles/2;

    for(int j=0; j<nParticles/2; j++) { //Cols
        for(int l=0; l<nParticles/2; l++) { //Rows
            sumSjUp(j) += newStates(l) * invSlaterMatrixUp(l,j);
            sumSjDown(j) += newStates(l) * invSlaterMatrixDown(l,j);
        }
    }

    //Update inverse matrices:

    //All columns except column corresponding to particle i:
    if(i<nParticles/2) { //If particle i has spin up
        for (int j=0; j<nParticles/2; j++) {
            for(int k=0; k<nParticles/2; k++) {
                if(j != particle) invSlaterMatrixUp(k,j) = invSlaterMatrixUp(k,j) - (sumSjUp(j)/ratio)*invSlaterMatrixUp(k,particle);
            }
        }
    }
    else {  //If particle i has spin down
        for (int j=0; j<nParticles/2; j++) { //Cols, inv matrix
            for(int k=0; k<nParticles/2; k++) { //Rows, inv matrix
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
    vec gradient = zeros<vec>(nDimensions,1);
    vec invMatrix = zeros<vec>(nParticles/2,1);

    for(int j=0; j<nParticles/2; j++) {
        if(i<nParticles/2) { invMatrix(j) = invSlaterMatrixUp(j,i); }
        else { invMatrix(j) = invSlaterMatrixDown(j,i-nParticles/2); }
    }

    for (int d=0; d<nDimensions; d++) { rtot += r(i,d)*r(i,d); }
    rtot = sqrt(rtot);

    gradient = function->dPsi1s(rtot, i, r, alpha)*(1/ratio)*invMatrix(0); //n=1,l=0,ml=0
    if(nParticles/2 > 1) gradient += function->dPsi2s(rtot, i, r, alpha)*(1/ratio)*invMatrix(1); //n=2,l=0,ml=0
    if(nParticles/2 > 2) gradient += function->dPsi2p_1(rtot, i, r, alpha)*(1/ratio)*invMatrix(2); //n=2,l=1,ml=-1
    if(nParticles/2 > 3) gradient += function->dPsi2p0(rtot, i, r, alpha)*(1/ratio)*invMatrix(3); //n=2,l=1,ml=0
    if(nParticles/2 > 4) gradient += function->dPsi2p1(rtot, i, r, alpha)*(1/ratio)*invMatrix(4); //n=2,l=1,ml=1


    //cout <<"Analytical: " << gradient<< endl;

    return gradient;

}

double slaterDeterminant::dWaveFunction_dalpha(const mat &r, double alpha) {

    double rtot = 0;
    double dalpha = 0;
    vec invMatrix = zeros<vec>(nParticles/2,1);

    for(int i=0; i<nParticles; i++) {

        invMatrix.fill(0);
        for(int j=0; j<nParticles/2; j++) {
            if(i<nParticles/2) { invMatrix(j) = invSlaterMatrixUp(j,i); }
            else { invMatrix(j) = invSlaterMatrixDown(j,i-nParticles/2); }
        }

        rtot = 0;
        for (int d=0; d<nDimensions; d++) { rtot += r(i,d)*r(i,d); }
        rtot = sqrt(rtot);

        dalpha += function->dPsi1s_dalpha(rtot, alpha)*invMatrix(0); //n=1,l=0,ml=0
        if(nParticles/2 > 1) dalpha += function->dPsi2s_dalpha(rtot, alpha)*invMatrix(1); //n=2,l=0,ml=0
        if(nParticles/2 > 2) dalpha += function->dPsi2p_1_dalpha(rtot, i, r, alpha)*invMatrix(2); //n=2,l=1,ml=-1
        if(nParticles/2 > 3) dalpha += function->dPsi2p0_dalpha(rtot, i, r, alpha)*invMatrix(3); //n=2,l=1,ml=0
        if(nParticles/2 > 4) dalpha += function->dPsi2p1_dalpha(rtot, i, r, alpha)*invMatrix(4); //n=2,l=1,ml=1
    }

    return dalpha;

}




vec slaterDeterminant::gradientWaveFunctionNum(const mat &r, int i, double alpha_, double beta_) {

vec grad = zeros(nDimensions);
mat rPlus = zeros<mat>(nDimensions);
mat rMinus = zeros<mat>(nDimensions);
double wf = beryllium(r,alpha_);
double h = 0.001;

rPlus = rMinus = r;

double waveFunctionMinus = 0;
double waveFunctionPlus = 0;

//First derivative

    for(int j = 0; j < nDimensions; j++) {
        rPlus(i,j) = r(i,j)+h;
        rMinus(i,j) = r(i,j)-h;

        waveFunctionMinus = beryllium(rMinus,alpha_);
        waveFunctionPlus = beryllium(rPlus,alpha_);
        grad(j) = (waveFunctionPlus - waveFunctionMinus)/(2*wf*h);
        rPlus(i,j) = r(i,j);
        rMinus(i,j) = r(i,j);
    }

//cout <<"Numerical: "<< grad << endl;

return grad;
}


double slaterDeterminant::laPlaceWaveFunction(const mat &r, double alpha, double beta) {

    double rtot = 0;
    double laplace = 0;

    vec invMatrix = zeros<vec>(nParticles/2,1);

    for(int i = 0; i < nParticles; i++) { //Particles

        rtot = 0;
        for (int d=0; d<nDimensions; d++) { rtot += r(i,d)*r(i,d); } //Get r for particle
        rtot = sqrt(rtot);

        invMatrix.fill(0);
        for(int j=0; j<nParticles/2; j++) {
            if(i<nParticles/2) { invMatrix(j) = invSlaterMatrixUp(j,i); }
            else { invMatrix(j) = invSlaterMatrixDown(j,i-nParticles/2); }
        }

        //Get laplacian for each state:
        laplace += function->d2Psi1s(rtot, alpha)*invMatrix(0); //n=1,l=0,ml=0
        if(nParticles/2 > 1) laplace += function->d2Psi2s(rtot, alpha)*invMatrix(1);//n=2,l=0,ml=0
        if(nParticles/2 > 2) laplace += function->d2Psi2p_1(rtot, i, r, alpha)*invMatrix(2);//n=2,l=1,ml=-1
        if(nParticles/2 > 3) laplace += function->d2Psi2p0(rtot, i, r, alpha)*invMatrix(3);//n=2,l=1,ml=0
        if(nParticles/2 > 4) laplace += function->d2Psi2p1(rtot, i, r, alpha)*invMatrix(4);//n=2,l=1,ml=1

    }

     //cout <<"Analytical: "<<laplace<<endl;

    return laplace;

}


double slaterDeterminant::laPlaceWaveFunctionNum(const mat &r, double alpha, double beta) {

    double h2 = 1000000;
    double h = 0.001;

    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    double waveFunctionCurrent = beryllium(r, alpha); //Find wavefunction for r

    //Second derivative (del^2):

    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = beryllium(rMinus, alpha);
            waveFunctionPlus = beryllium(rPlus, alpha);
            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
    kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

    //cout <<"Numerical: "<<kineticEnergy<<endl;

    return kineticEnergy;


}


double slaterDeterminant::beryllium(const mat &r, double &alpha_)  {

    double rs[nParticles];
    double sum = 0;

    //Find |r| for each electron:
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        rs[i] = sqrt(rSingleParticle);
        sum += rs[i];
    }

    //Slater determinant, Be
    double waveFunction = (function->psi1s(rs[0], alpha_)*function->psi2s(rs[1], alpha_) - function->psi1s(rs[1], alpha_)*function->psi2s(rs[0], alpha_))*
                          (function->psi1s(rs[2], alpha_)*function->psi2s(rs[3], alpha_) - function->psi1s(rs[3], alpha_)*function->psi2s(rs[2], alpha_));



    return waveFunction;

}


