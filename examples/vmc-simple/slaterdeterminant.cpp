#include "slaterdeterminant.h"

#include <armadillo>


using namespace arma;
using namespace std;

slaterDeterminant::slaterDeterminant(int nParticles_, int nDimensions_):

    nDimensions(nDimensions_),
    nParticles(nParticles_)

{
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
