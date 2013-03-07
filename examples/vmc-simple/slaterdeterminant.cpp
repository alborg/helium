#include "slaterdeterminant.h"

#include <armadillo>


using namespace arma;
using namespace std;

slaterDeterminant::slaterDeterminant(int nParticles_):

    nParticles(nParticles_)

{
}


double slaterDeterminant::beryllium(const mat &r, double &alpha_)  {

    double rs[nParticles];
    alpha = alpha_;

    for(int i=0; i<nParticles; i++) {
        rs[i] = norm(r.row(i), 2);
    }

    double waveFunction = (psi1s(rs[0])*psi2s(rs[1]) - psi1s(rs[1])*psi2s(rs[0]))*
                          (psi1s(rs[2])*psi2s(rs[3]) - psi1s(rs[3])*psi2s(rs[2]));

    return waveFunction;

}


double slaterDeterminant::psi1s(double r) {

    return exp(-alpha*r);

}



double slaterDeterminant::psi2s(double r) {

    return (1-(alpha*r)/2)*exp(-alpha*r/2);

}
