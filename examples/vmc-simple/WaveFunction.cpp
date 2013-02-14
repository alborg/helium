#include "WaveFunction.h"
#include "lib.h"

#include <armadillo>


using namespace arma;
using namespace std;

WaveFunction::WaveFunction(int nParticles_, int nDimensions_, double alpha_, double beta_) :

    nDimensions(nDimensions_),
    nParticles(nParticles_),
    alpha(alpha_),
    beta(beta_)


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
    double r12norm = norm(r12, 2);
    return exp(r12norm / (2 * (1 + beta * r12norm)));
}
