#include "WaveFunction.h"
#include "lib.h"

#include <armadillo>


using namespace arma;
using namespace std;

WaveFunction::WaveFunction(int nP, int nD) :

    nDimensions(nD),
    nParticles(nP)


{
}

double WaveFunction::waveFunction(const mat &r) {
    double argument = 0;
    for(int i = 0; i < nParticles; i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        argument += sqrt(rSingleParticle);
    }
    return exp(-argument * alpha) * jastrowFactor(r);
}


double WaveFunction::jastrowFactor(const mat &r) {
    rowvec r12 = r.row(1) - r.row(0);
    r12norm = norm(r12, 2);
    return exp(r12norm / (2 * (1 + beta * r12norm)));
}
