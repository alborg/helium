#include <unittest++/UnitTest++.h>
#include "../../vmc-simple/WaveFunction.h"
#include "../../vmc-simple/slaterdeterminant.h"
#include "../../vmc-simple/correlation.h"
#include "../../vmc-simple/lib.h"
#include <armadillo>

using namespace arma;
using namespace std;


TEST(Wave) {
    int nParticles = 4;
    int nDimensions = 3;
    double alpha = 4;
    double beta = 0.2;
    int stepLength = 1;
    long idum = -1;
    mat rOld = zeros<mat>(nParticles, nDimensions);
    mat rNew = zeros<mat>(nParticles,nDimensions);

    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
        }
    }

    rNew = rOld;

    WaveFunction *function = new WaveFunction(nParticles, nDimensions);
    slaterDeterminant *slater = new slaterDeterminant(nParticles, nDimensions);
    correlation *corr = new correlation(nParticles,nDimensions);

    slater -> buildDeterminant(rOld,alpha,beta);
    double test = slater->getDeterminant();

//    CHECK(slater->getDeterminant() == slater->beryllium(rOld,alpha));
//    CHECK(slater->getInvDeterminant() == slater->beryllium(rOld,alpha));

    for(int j=0; j<1; j++) { //Cycles
        for(int i = 0; i < nParticles; i++) { //Particles

            cout<<"cycle, particle: "<<j<<" "<<i<<endl;
            vec gradCorr = corr->gradientWaveFunction(rOld,i,beta);
            vec gradCorrNum = corr->gradientWaveFunctionNum(rOld,i,beta);
            cout <<"Grad Jastrow: "<< gradCorr<<" "<<gradCorrNum <<endl;
            double g_x = abs(gradCorr(0) - gradCorrNum(0));
            double g_y = abs(gradCorr(1) - gradCorrNum(1));
            double g_z = abs(gradCorr(2) - gradCorrNum(2));
            CHECK(g_x < 1e-2);
            //CHECK(g_y < 1e-2);
            //CHECK(g_z < 1e-2);

            double laplaceNum = corr->laPlaceWaveFunctionNum(rOld, beta);
            double laplace = corr->laPlaceWaveFunction(rOld, beta);
            double deltaLaP = abs(laplace - laplaceNum);
            CHECK(laplaceNum < 1e-2);

            //        vec a1 = slater->gradientWaveFunction(rOld,i,1.0,4,0.2);
            //        vec b1 = slater->gradientWaveFunctionNum(rOld,i,4,0.2);

            //        double error0 = 1e-2;
            //        double gradient_x1 = abs(a1(0) - b1(0));
            //        CHECK(gradient_x1 < error0);


            //New position:
            for(int d = 0; d < nDimensions; d++) {
                rNew(i,d) = rOld(i,d) + stepLength*(ran2(&idum) - 0.5);
            }

            //Accept new step
            //        double ratio1 = slater->getRatioDeterminantNum(i,rOld,rNew,alpha, beta);
            //        double ratio2 = slater->getRatioDeterminant(i,rNew,alpha, beta);

            //        double error3 = 1e-3;
            //        double getRatioDeterminant = abs(ratio2 - ratio1);
            //        CHECK(getRatioDeterminant < error3);

            //        vec a = slater->gradientWaveFunction(rNew,i,ratio1,4,0.2);
            //        vec b = slater->gradientWaveFunctionNum(rNew,i,4,0.2);

            double error2 = 1e-2;
            //        double gradient_x = abs(a(0) - b(0));
            //        CHECK(gradient_x < error2);

            //        double gradient_y = abs(a(1) - b(1));
            //        double gradient_z = abs(a(2) - b(2));
            //        CHECK(gradient_y < error2);
            //        CHECK(gradient_z < error2);



            //        slater->updateDeterminant(rNew, rOld, i, alpha, beta, ratio1);

            //        double error4 = 1e-3;
            //        double getDeterminant = abs(slater->getInvDeterminant() - slater->beryllium(rNew,alpha));

            //        CHECK(getDeterminant < error4);

                    for(int d = 0; d < nDimensions; d++) {
                        rOld(i,d) = rNew(i,d);
                    }

            //        double lpN = slater->laPlaceWaveFunctionNum(rOld,alpha,beta);
            //        double lp = slater->laPlaceWaveFunction(rOld,alpha,beta);
            //        double lap_diff = abs(lpN - lp);
            //        double error1 = 1e-2;
            //        //cout << lpN/lp<<endl;
            //        CHECK(lap_diff<error1);
        }
    }
}


int main()
{
    return UnitTest::RunAllTests();
}
