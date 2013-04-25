#include <unittest++/UnitTest++.h>
#include "../../vmc-simple/WaveFunction.h"
#include "../../vmc-simple/slaterdeterminant.h"
#include "../../vmc-simple/lib.h"
#include <armadillo>

using namespace arma;

TEST(Wave) {
    int nParticles = 4;
    int nDimensions = 3;
    double alpha = 4;
    double beta = 0.2;
    int stepLength = 1;
    long idum = -1.0;
    mat rOld = zeros<mat>(nParticles, nDimensions);
    mat rNew = zeros<mat>(nParticles,nDimensions);

    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
        }
    }

    rNew = rOld;

    slaterDeterminant *slater = new slaterDeterminant(nParticles, nDimensions);
    slater -> buildDeterminant(rOld,alpha,beta);

    CHECK(slater->getDeterminant() == slater->beryllium(rOld,alpha));

    for(int i = 0; i < nParticles; i++) { //Particles

        //New position:
        for(int d = 0; d < nDimensions; d++) {
            rNew(i,d) = rOld(i,d) + stepLength*(ran2(&idum) - 0.5);
        }


        vec a = zeros<vec>(nDimensions,1);
        a = slater->gradientWaveFunction(rOld,i,1,4,0.2);
        vec b  = zeros<vec>(nDimensions,1);
        b = slater->gradientWaveFunctionNum(rOld,i,4,0.2);

//        double gradient_x = abs(a(0) - b(0));
//        double gradient_y = abs(a(1) - b(1));
//        double gradient_z = abs(a(2) - b(2));
//        double error2 = 1e-3;
//        CHECK(gradient_x < error2);
//        CHECK(gradient_y < error2);
//        CHECK(gradient_z < error2);

        double error3 = 1e-3;
        double getRatioDeterminant = abs(slater->getRatioDeterminant(i,rOld,alpha, beta) - slater->getRatioDeterminantNum(i,rOld,rNew,alpha, beta));
        CHECK(getRatioDeterminant < error3);

        //Accept new step
        //double ratio = slater->getRatioDeterminantNum(i,rOld,rNew,alpha, beta);
        double ratio = slater->getRatioDeterminant(i,rOld,alpha, beta);
        slater->updateDeterminant(rNew, rOld, i, alpha, beta, ratio);

        double error4 = 1e-3;
        double getDeterminant = abs(slater->getInvDeterminant() - slater->beryllium(rNew,alpha));
        CHECK(getDeterminant < error4);

        for(int d = 0; d < nDimensions; d++) {
            rOld(i,d) = rNew(i,d);
        }
    }
}


int main()
{
    return UnitTest::RunAllTests();
}
