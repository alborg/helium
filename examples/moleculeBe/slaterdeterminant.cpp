#include "slaterdeterminant.h"
#include "WaveFunction.h"

#include <armadillo>


using namespace arma;
using namespace std;

slaterDeterminant::slaterDeterminant(int nDimensions_, int nProtons_, int nElectrons_):

    nDimensions(nDimensions_),
    nProtons(nProtons_),
    nElectrons(nElectrons_),
    nParticles(nElectrons*nProtons),
    slaterMatrixUp1(zeros(nElectrons/2,nElectrons/2)),
    slaterMatrixDown1(zeros(nElectrons/2,nElectrons/2)),
    slaterMatrixUp2(zeros(nElectrons/2,nElectrons/2)),
    slaterMatrixDown2(zeros(nElectrons/2,nElectrons/2)),
    invSlaterMatrixUp1(zeros(nElectrons/2,nElectrons/2)),
    invSlaterMatrixDown1(zeros(nElectrons/2,nElectrons/2)),
    invSlaterMatrixUp2(zeros(nElectrons/2,nElectrons/2)),
    invSlaterMatrixDown2(zeros(nElectrons/2,nElectrons/2)),
    function(new WaveFunction(nParticles,nDimensions))

{
}


void slaterDeterminant::buildDeterminant(const mat &r, double &alpha_) {

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

    //State 1s
    slaterMatrixUp1(0,0) = function->psi1s(rs[0], alpha);
    slaterMatrixUp1(1,0) = function->psi1s(rs[1], alpha);

    slaterMatrixDown1(0,0) = function->psi1s(rs[2], alpha);
    slaterMatrixDown1(1,0) = function->psi1s(rs[3], alpha);

    slaterMatrixUp2(0,0) = function->psi1s(rs[4], alpha);
    slaterMatrixUp2(1,0) = function->psi1s(rs[5], alpha);

    slaterMatrixDown2(0,0) = function->psi1s(rs[6], alpha);
    slaterMatrixDown2(1,0) = function->psi1s(rs[7], alpha);

    //State 2s
    slaterMatrixUp1(0,1) = function->psi2s(rs[0], alpha);
    slaterMatrixUp1(1,1) = function->psi2s(rs[1], alpha);

    slaterMatrixDown1(0,1) = function->psi2s(rs[2], alpha);
    slaterMatrixDown1(1,1) = function->psi2s(rs[3], alpha);

    slaterMatrixUp2(0,1) = function->psi2s(rs[4], alpha);
    slaterMatrixUp2(1,1) = function->psi2s(rs[5], alpha);

    slaterMatrixDown2(0,1) = function->psi2s(rs[6], alpha);
    slaterMatrixDown2(1,1) = function->psi2s(rs[7], alpha);

cout<<rs<<endl;

    invSlaterMatrixUp1 = inv(slaterMatrixUp1);
    invSlaterMatrixDown1 = inv(slaterMatrixDown1);
    invSlaterMatrixUp2 = inv(slaterMatrixUp2);
    invSlaterMatrixDown2 = inv(slaterMatrixDown2);

}


double slaterDeterminant::getDeterminant() {

    return det(slaterMatrixUp1)*det(slaterMatrixDown1)*det(slaterMatrixUp2)*det(slaterMatrixDown2);

}


double slaterDeterminant::getInvDeterminant() {

    return 1/(det(invSlaterMatrixUp1)*det(invSlaterMatrixDown1)) * 1/(det(invSlaterMatrixUp2)*det(invSlaterMatrixDown2));

}


double slaterDeterminant::getRatioDeterminant(int i, const mat &r, double alpha) {

    double ratio = 0;
    double rSingleParticle = 0;
    //Get rtot for particle i's new position:
    for(int d = 0; d < nDimensions; d++) {
        rSingleParticle += r(i,d) * r(i,d);
    }
    double rtot = sqrt(rSingleParticle);

    vec updatedStates = getStates(rtot, alpha); //Get states with new r for particle i

    if(i==0 || i==1) { //Atom 1, particle spin up
        for(int j=0; j<nElectrons/2; j++) { //States
            ratio += updatedStates(j) * invSlaterMatrixUp1(j,i);
        }
    }
    if(i==2 || i==3) { //Atom 1, particle spin down
        for(int j=0; j<nElectrons/2; j++) { //States
            ratio += updatedStates(j) * invSlaterMatrixDown1(j,i-2);
        }
    }
    if(i==4 || i==5) { //Atom 2, particle spin up
        for(int j=0; j<nElectrons/2; j++) { //States
            ratio += updatedStates(j) * invSlaterMatrixUp2(j,i-4);
        }
    }
    if(i==6 || i==7) { //Atom 2, particle spin down
        for(int j=0; j<nElectrons/2; j++) { //States
            ratio += updatedStates(j) * invSlaterMatrixDown2(j,i-6);
        }
    }



    return ratio;
}


double slaterDeterminant::getRatioDeterminantNum(int i, const mat &rOld, const mat &rNew, double alpha) {

    return beryllium(rNew, alpha) / beryllium(rOld, alpha);

}

vec slaterDeterminant::getStates(double rtot, double alpha) {

    vec updatedStates = zeros<vec>(nElectrons/2,1);

    updatedStates(0) = function->psi1s(rtot, alpha); //n=1,l=0,ml=0
    updatedStates(1) = function->psi2s(rtot, alpha); //n=2,l=0,ml=0

    return updatedStates;
}



void slaterDeterminant::updateDeterminant(const mat &rNew, const mat &rOld, int i, double &alpha_, double ratio) {

    double rtot = 0;
    int particle = 0;
    //Get rtot for particles' position:
    for(int d = 0; d < nDimensions; d++) {
        rtot += rNew(i,d) * rNew(i,d);
    }
    rtot = sqrt(rtot);

    vec newStates = getStates(rtot, alpha_);

    vec sumSj = zeros<vec>(nElectrons/2); //Sum over states(l)*d_lj for particles j

    if(i<nElectrons) { //Atom1, electrons 0-3

        particle = i;
        if(i>1) particle = i - nElectrons/2;

        if(i==0 || i==1) { //Spin up

            for(int j=0; j<nElectrons/2; j++) { //Cols
                for(int l=0; l<nElectrons/2; l++) { //Rows
                    sumSj(j) += newStates(l) * invSlaterMatrixUp1(l,j);
                }
            }

            //Update all columns except column corresponding to particle i:

            for (int j=0; j<nElectrons/2; j++) {
                for(int k=0; k<nElectrons/2; k++) {
                    if(j != i) invSlaterMatrixUp1(k,j) = invSlaterMatrixUp1(k,j) - (sumSj(j)/ratio)*invSlaterMatrixUp1(k,i);
                }
            }
        }
        else {  //If particle i has spin down

            for(int j=0; j<nElectrons/2; j++) { //Cols
                for(int l=0; l<nElectrons/2; l++) { //Rows
                    sumSj(j) += newStates(l) * invSlaterMatrixDown1(l,j);
                }
            }

            for (int j=0; j<nElectrons/2; j++) { //Cols, inv matrix
                for(int k=0; k<nElectrons/2; k++) { //Rows, inv matrix
                    if(j != particle) invSlaterMatrixDown1(k,j) = invSlaterMatrixDown1(k,j) - (sumSj(j)/ratio)*invSlaterMatrixDown1(k,particle);
                }
            }
        }

        //Update column corresponding to particle i:
        for(int k=0; k<nElectrons/2; k++) { //States (rows)
            if(i<nElectrons/2) { invSlaterMatrixUp1(k,particle) = (1/ratio)*invSlaterMatrixUp1(k,particle); }
            else {invSlaterMatrixDown1(k,particle) = (1/ratio)*invSlaterMatrixDown1(k,particle); }
        }

    }
    else { //Atom2, electrons 4-7

        i -= 4;
        particle = i;
        if(i>1) particle = i - nElectrons/2;

        if(i==0 || i==1) { //Spin up

            for(int j=0; j<nElectrons/2; j++) { //Cols
                for(int l=0; l<nElectrons/2; l++) { //Rows
                    sumSj(j) += newStates(l) * invSlaterMatrixUp2(l,j);
                }
            }

            //Update all columns except column corresponding to particle i:

            for (int j=0; j<nElectrons/2; j++) {
                for(int k=0; k<nElectrons/2; k++) {
                    if(j != i) invSlaterMatrixUp2(k,j) = invSlaterMatrixUp2(k,j) - (sumSj(j)/ratio)*invSlaterMatrixUp2(k,i);
                }
            }
        }
        else {  //If particle i has spin down

            for(int j=0; j<nElectrons/2; j++) { //Cols
                for(int l=0; l<nElectrons/2; l++) { //Rows
                    sumSj(j) += newStates(l) * invSlaterMatrixDown2(l,j);
                }
            }

            for (int j=0; j<nElectrons/2; j++) { //Cols, inv matrix
                for(int k=0; k<nElectrons/2; k++) { //Rows, inv matrix
                    if(j != particle) invSlaterMatrixDown2(k,j) = invSlaterMatrixDown2(k,j) - (sumSj(j)/ratio)*invSlaterMatrixDown2(k,particle);
                }
            }
        }

        //Update column corresponding to particle i:
        for(int k=0; k<nElectrons/2; k++) { //States (rows)
            if(i<nElectrons/2) { invSlaterMatrixUp2(k,particle) = (1/ratio)*invSlaterMatrixUp2(k,particle); }
            else {invSlaterMatrixDown2(k,particle) = (1/ratio)*invSlaterMatrixDown2(k,particle); }
        }
    }



}

vec slaterDeterminant::gradientWaveFunction(const mat &r, int i, double ratio, double alpha) {

    vec gradient = zeros<vec>(nDimensions,1);
    vec invMatrix = zeros<vec>(nElectrons/2,1);

    for(int j=0; j<nElectrons/2; j++) {
        if(i==0 || i==1) { invMatrix(j) = invSlaterMatrixUp1(j,i); }
        if(i==2 || i==3) { invMatrix(j) = invSlaterMatrixDown1(j,i-2); }
        if(i==4 || i==5) { invMatrix(j) = invSlaterMatrixUp2(j,i-4); }
        if(i==6 || i==7) { invMatrix(j) = invSlaterMatrixDown2(j,i-6); }
    }

    double rtot = 0;
    for (int d=0; d<nDimensions; d++) { rtot += r(i,d)*r(i,d); }
    rtot = sqrt(rtot);
    cout<<i<<" "<<rtot<<" "<<invMatrix<<endl;

    gradient = function->dPsi1s(rtot, i, r, alpha)*(1/ratio)*invMatrix(0); //n=1,l=0,ml=0
    gradient += function->dPsi2s(rtot, i, r, alpha)*(1/ratio)*invMatrix(1); //n=2,l=0,ml=0

    //cout <<"Analytical: " << gradient<< endl;

    return gradient;

}




vec slaterDeterminant::gradientWaveFunctionNum(const mat &r, int i, double alpha_) {

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


double slaterDeterminant::laPlaceWaveFunction(const mat &r, int i, double alpha) {

    double rtot = 0;
    double laplace = 0;

    vec invMatrix = zeros<vec>(nElectrons/2,1);


    rtot = 0;
    for (int d=0; d<nDimensions; d++) { rtot += r(i,d)*r(i,d); } //Get r for particle
    rtot = sqrt(rtot);

    invMatrix.fill(0);
    for(int j=0; j<nElectrons/2; j++) {
        if(i==0 || i==1) { invMatrix(j) = invSlaterMatrixUp1(j,i); }
        if(i==2 || i==3) { invMatrix(j) = invSlaterMatrixDown1(j,i-2); }
        if(i==4 || i==5) { invMatrix(j) = invSlaterMatrixUp2(j,i-4); }
        if(i==6 || i==7) { invMatrix(j) = invSlaterMatrixDown2(j,i-6); }
    }

    //Get laplacian for each state:
    laplace += function->d2Psi1s(rtot, alpha)*invMatrix(0); //n=1,l=0,ml=0
    laplace += function->d2Psi2s(rtot, alpha)*invMatrix(1);//n=2,l=0,ml=0



    //cout <<"Analytical: "<<laplace<<endl;

    return laplace;

}


double slaterDeterminant::laPlaceWaveFunctionNum(const mat &r, double alpha) {

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
    double waveFunction1 = (function->psi1s(rs[0], alpha_)*function->psi2s(rs[1], alpha_) - function->psi1s(rs[1], alpha_)*function->psi2s(rs[0], alpha_))*
                          (function->psi1s(rs[2], alpha_)*function->psi2s(rs[3], alpha_) - function->psi1s(rs[3], alpha_)*function->psi2s(rs[2], alpha_));

    double waveFunction2 = (function->psi1s(rs[4], alpha_)*function->psi2s(rs[5], alpha_) - function->psi1s(rs[5], alpha_)*function->psi2s(rs[4], alpha_))*
                          (function->psi1s(rs[6], alpha_)*function->psi2s(rs[7], alpha_) - function->psi1s(rs[7], alpha_)*function->psi2s(rs[6], alpha_));


    return waveFunction1*waveFunction2;

}


