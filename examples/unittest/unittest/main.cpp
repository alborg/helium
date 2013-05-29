#include <unittest++/UnitTest++.h>
//#include "../../molecules/WaveFunction.h"
//#include "../../molecules/lib.h"
//#include "../../molecules/hamiltonian.h"
//#include "../../vmc-simple/vmcsolver.h"
//#include "../../vmc-simple/slaterdeterminant.h"
//#include "../../vmc-simple/correlation.h"
//#include "../../vmc-simple/hamiltonian.h"
#include "../../moleculeBe/vmcsolver.h"
#include "../../moleculeBe/slaterdeterminant.h"
#include "../../moleculeBe/correlation.h"
#include "../../moleculeBe/hamiltonian.h"
#include "../../moleculeBe/lib.h"

#include <armadillo>

using namespace arma;
using namespace std;


TEST(Wave) {
    int nParticles = 2;
    int charge = 8;
    int nDimensions = 3;
    double alpha = 3.7;
    double beta = 0.23;
    int stepLength = 1;
    long idum = -1;
    double h = 0.001;
    double h2 = 1000000;

    double nProtons = 2;
    double nElectrons = 4;
    double R = 4.63;

    nParticles = 8;
    mat rOld = zeros<mat>(nParticles, nDimensions);
    mat rNew = zeros<mat>(nParticles,nDimensions);

    mat rProtons = zeros<mat>(nProtons, nDimensions);
    rProtons(0,2) = -R/2;
    rProtons(1,2) = R/2;

    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
        }
    }

    rNew = rOld;

    WaveFunction *function = new WaveFunction(nParticles,nDimensions);
    slaterDeterminant *slater = new slaterDeterminant(nDimensions, nProtons, nElectrons);
    correlation *corr = new correlation(nDimensions, nProtons, nElectrons);
    Hamiltonian *hamiltonian = new Hamiltonian(nDimensions,h,h2,charge, nProtons, nElectrons, R, rProtons);

    slater -> buildDeterminant(rOld,alpha);
   // double test = slater->getDeterminant();
    cout<<"Determinants: "<< slater->getDeterminant()<<" "<<slater->beryllium(rOld,alpha)<<endl;

    //CHECK(slater->getDeterminant() == slater->beryllium(rOld,alpha));
////    CHECK(slater->getInvDeterminant() == slater->beryllium(rOld,alpha));

//    double EnSum=0;
//    double EHeSum=0;
    int nCycles = 5;

    for(int j=0; j<nCycles; j++) { //Cycles
        for(int i = 0; i < nParticles; i++) { //Particles

//            double En = hamiltonian->localEnergy(rOld,alpha,beta,slater,corr);
//            EnSum +=En;

//            double EHe = hamiltonian->analyticEnergyHe(rOld,alpha,beta);
//            EHeSum += EHe;



////           cout<<"En, EHe: "<<En<<" "<<EHe<<endl;


            cout<<"cycle, particle: "<<j<<" "<<i<<endl;

            //cout<<"Jastrow: "<<corr->jastrow(rOld,beta)<<" "<<corr->jastrowNum(rOld,beta)<<endl;

            vec gradCorr = corr->gradientWaveFunction(rOld,i,beta);
            vec gradCorrNum = corr->gradientWaveFunctionNum(rOld,i,beta);
//           // cout <<"Grad Jastrow: "<< gradCorr<<" "<<gradCorrNum <<endl;
            double g_x = abs(gradCorr(0) - gradCorrNum(0));
            double g_y = abs(gradCorr(1) - gradCorrNum(1));
            double g_z = abs(gradCorr(2) - gradCorrNum(2));
            CHECK(g_x < 1e-2);
            CHECK(g_y < 1e-2);
            CHECK(g_z < 1e-2);

            double laplaceNum = corr->laPlaceWaveFunctionNum(rOld, beta);
            double laplace = corr->laPlaceWaveFunction(rOld, beta);
            double deltaLaP = abs(laplace - laplaceNum);
            CHECK(laplaceNum < 1e-2);

            vec a1 = slater->gradientWaveFunction(rOld,i,1,alpha);
            vec b1 = slater->gradientWaveFunctionNum(rOld,i,alpha);
            cout<<"Gradient, slater: "<<a1<<" "<<b1<<endl;
            double error0 = 1e-2;
            double gradient_x1 = abs(a1(0) - b1(0));
            CHECK(gradient_x1 < error0);

            double lpN = slater->laPlaceWaveFunctionNum(rOld,alpha);
            double lp = slater->laPlaceWaveFunction(rOld,alpha);
            double lap_diff = abs(lpN + 0.5*lp);
            double error1 = 1e-2;
            cout <<"Laplace, slater: "<< -0.5*lp<<" "<<lpN<<endl;
            CHECK(lap_diff<error1);


            //New position:
            for(int d = 0; d < nDimensions; d++) {
                rNew(i,d) = rOld(i,d) + stepLength*(ran2(&idum) - 0.5);
            }

            //            //Accept new step
            double ratio1 = slater->getRatioDeterminantNum(i,rOld,rNew,alpha);
            double ratio2 = slater->getRatioDeterminant(i,rNew,alpha);

            cout<<"Ratio: "<<ratio2<<" "<<ratio1<<endl;

            double error3 = 1e-3;
            double getRatioDeterminant = abs(ratio2 - ratio1);
            CHECK(getRatioDeterminant < error3);

//            //        vec a = slater->gradientWaveFunction(rNew,i,ratio1,4,0.2);
//            //        vec b = slater->gradientWaveFunctionNum(rNew,i,4,0.2);

//            double error2 = 1e-2;
//            //        double gradient_x = abs(a(0) - b(0));
//            //        CHECK(gradient_x < error2);

//            //        double gradient_y = abs(a(1) - b(1));
//            //        double gradient_z = abs(a(2) - b(2));
//            //        CHECK(gradient_y < error2);
//            //        CHECK(gradient_z < error2);



                    slater->updateDeterminant(rNew, rOld, i, alpha, ratio1);

                    double error4 = 1e-3;
                    double test1 = slater->beryllium(rNew,alpha);
                    cout<<"New r: "<<slater->getInvDeterminant()<<" "<<test1<<endl;
                    double getDeterminant = abs(slater->getInvDeterminant() - test1);

                    CHECK(getDeterminant < error4);

                    for(int d = 0; d < nDimensions; d++) {
                        rOld(i,d) = rNew(i,d);
                    }


        }
    }

//        EnSum = EnSum/(nCycles*nParticles);
//        EHeSum = EHeSum/(nCycles*nParticles);
//        cout << "Analytic, num: "<<EHeSum<<" "<<EnSum<<endl;




}

//TEST(molecule) {

//    int stepLength = 1;
//    long idum = -1;

//    int nDimensions = 2;
//    int nProtons = 2;
//    int nElectrons = 4;
//    int nParticles = nElectrons*nProtons;
//    double R = 4.63;

//    double alpha = 3.5;
//    double beta = 0.3;

//    mat r = zeros<mat>(nParticles,nDimensions);
//    for(int i = 0; i < nParticles; i++) {
//        for(int j = 0; j < nDimensions; j++) {
//            r(i,j) = stepLength * (ran2(&idum) - 0.5);
//        }
//    }

//    mat rProtons = zeros<mat>(nProtons,nDimensions);
//    rProtons(0,2) = -R/2;
//    rProtons(1,2) = R/2;

//    WaveFunction *function = new WaveFunction(nDimensions,nProtons,nElectrons);
//    double wf = function->waveFunction(r,rProtons,alpha,beta);

//    rowvec r1p1 = r.row(0) - rProtons.row(0);
//    rowvec r1p2 = r.row(0) - rProtons.row(1);
//    rowvec r2p1 = r.row(1) - rProtons.row(0);
//    rowvec r2p2 = r.row(1) - rProtons.row(1);
//    rowvec r3p1 = r.row(2) - rProtons.row(0);
//    rowvec r3p2 = r.row(2) - rProtons.row(1);
//    rowvec r4p1 = r.row(3) - rProtons.row(0);
//    rowvec r4p2 = r.row(3) - rProtons.row(1);

//    rowvec r5p1 = r.row(4) - rProtons.row(0);
//    rowvec r5p2 = r.row(4) - rProtons.row(1);
//    rowvec r6p1 = r.row(5) - rProtons.row(0);
//    rowvec r6p2 = r.row(5) - rProtons.row(1);
//    rowvec r7p1 = r.row(6) - rProtons.row(0);
//    rowvec r7p2 = r.row(6) - rProtons.row(1);
//    rowvec r8p1 = r.row(7) - rProtons.row(0);
//    rowvec r8p2 = r.row(7) - rProtons.row(1);

//    double R1p1 = sqrt(pow(r1p1(0),2) + pow(r1p1(1),2) + pow(r1p1(2),2));
//    double R1p2 = sqrt(pow(r1p2(0),2) + pow(r1p2(1),2) + pow(r1p2(2),2));
//    double R2p1 = sqrt(pow(r2p1(0),2) + pow(r2p1(1),2) + pow(r2p1(2),2));
//    double R2p2 = sqrt(pow(r2p2(0),2) + pow(r2p2(1),2) + pow(r2p2(2),2));
//    double R3p1 = sqrt(pow(r3p1(0),2) + pow(r3p1(1),2) + pow(r3p1(2),2));
//    double R3p2 = sqrt(pow(r3p2(0),2) + pow(r3p2(1),2) + pow(r3p2(2),2));
//    double R4p1 = sqrt(pow(r4p1(0),2) + pow(r4p1(1),2) + pow(r4p1(2),2));
//    double R4p2 = sqrt(pow(r4p2(0),2) + pow(r4p2(1),2) + pow(r4p2(2),2));

//    double R5p1 = sqrt(pow(r5p1(0),2) + pow(r5p1(1),2) + pow(r5p1(2),2));
//    double R5p2 = sqrt(pow(r5p2(0),2) + pow(r5p2(1),2) + pow(r5p2(2),2));
//    double R6p1 = sqrt(pow(r6p1(0),2) + pow(r6p1(1),2) + pow(r6p1(2),2));
//    double R6p2 = sqrt(pow(r6p2(0),2) + pow(r6p2(1),2) + pow(r6p2(2),2));
//    double R7p1 = sqrt(pow(r7p1(0),2) + pow(r7p1(1),2) + pow(r7p1(2),2));
//    double R7p2 = sqrt(pow(r7p2(0),2) + pow(r7p2(1),2) + pow(r7p2(2),2));
//    double R8p1 = sqrt(pow(r8p1(0),2) + pow(r8p1(1),2) + pow(r8p1(2),2));
//    double R8p2 = sqrt(pow(r8p2(0),2) + pow(r8p2(1),2) + pow(r8p2(2),2));


//    double wfnum = (exp(-alpha*R1p1) + exp(-alpha*R1p2)) * (exp(-alpha*R2p1) + exp(-alpha*R2p2)) * (exp(-alpha*R3p1) + exp(-alpha*R3p2)) * (exp(-alpha*R4p1) + exp(-alpha*R4p2)) * (exp(-alpha*R5p1) + exp(-alpha*R5p2)) * (exp(-alpha*R6p1) + exp(-alpha*R6p2)) * (exp(-alpha*R7p1) + exp(-alpha*R7p2)) * (exp(-alpha*R8p1) + exp(-alpha*R8p2));
//    cout<<"num wave: "<<wfnum<<endl;

//    double jastrow = 0;
//    double a = 0.5;
//    rowvec r12;
//    double R12 = 0;
//    for(int i=1; i<nParticles; i++) {
//        for(int j=0; j<i; j++) {
//            r12 = r.row(i) - r.row(j);
//            R12 = sqrt(pow(r12(0),2) + pow(r12(1),2) + pow(r12(2),2));
//            if((i==1 && j==0) || (i==3 && j==2) || (i==4 && (j==0 || j==1)) || (i==5 && (j==0 || j==1 || j==4)) || (i==6 && (j==2 || j==3)) || (i==7 && (j==2 || j==3 || j==6))) {
//                a=0.25;
//            }
//            else a=0.5;
//            jastrow += a*R12/(1+beta*R12)/2;
//        }
//    }


//    jastrow = exp(jastrow);
//    wfnum = wfnum * jastrow;
//    cout <<"num jastrow: "<<jastrow<<endl;
//    double delta = abs(wfnum - wf);

//    cout <<"wavwfunc, num: "<<wf<<" "<<wfnum<<endl;

//    CHECK(delta < 1e-3);
//}

//TEST(dEdP) {

//    int nParticles = 4;
//    int charge = 4;
//    int nDimensions = 3;
//    long idum = -1;
//    double alpha_min = 3;
//    double alpha_max = 4;
//    double alpha_steps = 11;
//    double beta_min = 0.2;
//    double beta_max = 0.3;
//    double beta_steps = 11;
//    int steps = 10000;
//    double h = 0.001;
//    double h2 = 1000000;
//    double* allEnergies = new double[steps+1];
//    vec dE = zeros<vec>(2);
//    vec dPsi = zeros<vec>(2,1);
//    vec dPsi_Elocal = zeros<vec>(2,1);
//    double da = 0.2;
//    double db = 0.05;

//     slaterDeterminant *slater = new slaterDeterminant(nParticles, nDimensions);
//     correlation *corr = new correlation(nParticles,nDimensions);
//     Hamiltonian *hamiltonian = new Hamiltonian(nParticles,nDimensions,h,h2,charge);
//     VMCSolver *solver = new VMCSolver();

//     mat Emat = zeros<mat>(alpha_steps,beta_steps);
//     mat dEmat_alpha = zeros<mat>(alpha_steps,beta_steps);
//     mat dEmat_beta = zeros<mat>(alpha_steps,beta_steps);
//     vec alphas = zeros<vec>(alpha_steps);
//     vec betas = zeros<vec>(beta_steps);

//     double alpha = 0;
//     double beta = 0;

//     double alpha_step = (alpha_max - alpha_min)/(alpha_steps-1);
//     double beta_step = (beta_max - beta_min)/(beta_steps-1);
//     if(alpha_max == alpha_min)  alpha_step = 1;
//     if(beta_max == beta_min) beta_step = 1;


//     for(int i=0; i<alpha_steps; i++) {
//         alpha =  alpha_min + i*alpha_step;
//         alphas(i) = alpha;
//         for(int j=0; j<beta_steps; j++) {
//             beta = beta_min + j*beta_step;
//             betas(j) = beta;
//             vec Es = solver->MCImportance(idum, alpha,beta,steps, slater, hamiltonian, corr, allEnergies);
//             double E = Es(0)/(steps * nParticles);
//             Emat(i,j) = E;
//             dPsi(0) = Es(2)/(steps * nParticles);
//             dPsi(1) = Es(3)/(steps * nParticles);
//             dPsi_Elocal(0) = Es(4)/(steps * nParticles);
//             dPsi_Elocal(1) = Es(5)/(steps * nParticles);
//             dE = solver->gradE(dPsi, E, dPsi_Elocal);
//             dEmat_alpha(i,j) = dE(0);
//             dEmat_beta(i,j) = dE(1);
//             cout <<"Alpha, beta, E, dE: "<<alpha<<" "<< beta<<" "<<E<<" "<<dE(0)<<" "<<dE(1) <<endl;
//         }
//     }

//     cout <<alphas<<endl;
//     cout<<betas<<endl;

//     cout <<Emat<<endl;
//     cout <<dEmat_alpha + dEmat_beta<<endl;
//    cout <<dEmat_alpha<<endl;
//    cout <<dEmat_beta<<endl;


//}



int main()
{
    return UnitTest::RunAllTests();
}
