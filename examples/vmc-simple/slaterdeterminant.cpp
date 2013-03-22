#include "slaterdeterminant.h"

#include <armadillo>


using namespace arma;
using namespace std;

slaterDeterminant::slaterDeterminant(int nParticles_, int nDimensions_):

    nDimensions(nDimensions_),
    nParticles(nParticles_)

{
}


double slaterDeterminant::determinant(const mat &r, double &alpha_) {

    double rs[nParticles];
    alpha = alpha_; //Set class variable alpha

    mat spinUp = zeros(nParticles/2,nParticles/2);
    mat spinDown = zeros(nParticles/2,nParticles/2);

    int nElectrons[3] = {1,1,3}; //Half of 2,2,6, ignoring spin

    //Find |r| for each electron:
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        rs[i] = sqrt(rSingleParticle);
    }


    bool full = false;
    int state = 0;
    int totalStates = 0;
    int totalElectrons = 0;
    int endLoop = 0;
    double up = 0;
    double down = 0;

    do{
        if(nParticles/2 >= (totalElectrons + nElectrons[state])) {
            if (nParticles/2 == (totalElectrons + nElectrons[state])) full = true;
            totalElectrons += nElectrons[state];
            endLoop = nElectrons[state];

        }
        else {
            full = true;
            int remainderE = nParticles/2 - totalElectrons;
            endLoop = remainderE;
        }

        for(int i = 0; i < nParticles/2; i++) {

            if (state == 0) {up = psi1s(rs[i]); down = psi1s(rs[nParticles/2+i]);}
            if (state == 1) {up = psi2s(rs[i]); down = psi2s(rs[nParticles/2+i]);}
            if (state == 2) {up = psi2p(rs[i]); down = psi2p(rs[nParticles/2+i]);}

            for (int k = 0; k <endLoop ; k++) {
                spinUp(i,totalStates + k) = up;
                spinDown(i,totalStates + k) = down;
            }
        }
        totalStates += endLoop;

        ++state;

    }while(!full);


    double waveFunction = det(spinUp)*det(spinDown);
    return waveFunction;

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

    //return pow(alpha,6)*waveFunction/(4*sqrt(2)*acos(-1));
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
