#include "correlation.h"
#include <armadillo>

using namespace arma;
using namespace std;

correlation::correlation(int nParticles_, int nDimensions_):

    nDimensions(nDimensions_),
    nParticles(nParticles_)
{
}


double correlation::jastrow(const mat &r, double beta) {

    double exponent=0;
    double a = 0.25;
    double rij = 0;

    for(int j=1;j<nParticles;j++) {
        for(int i=0;i<j;i++) {
            int test1 = i-nParticles/2;
            int test2 = j-nParticles/2;
            if(test1*test2>0) a=0.25;
            else a=0.5;
            rij = 0;
            for(int d=0;d<nDimensions;d++) rij += (r(i,d)-r(j,d))*(r(i,d)-r(j,d));
            rij = sqrt(rij);
            exponent += (a*rij)/(1.0+beta*rij);
        }
    }

    return exp(exponent);

}

double correlation::jastrowNum(const mat &r, double beta) {

    double r12 = 0;
    double wf = 0;

    for(int i=0; i<nParticles; i++) {
        for(int p=0; p<i; p++) {
            r12 = 0;
            for(int d=0;d<nDimensions;d++) r12 += pow(r(i,d)-r(p,d),2);
            r12 = sqrt(r12);
            wf += r12/(2*(1+beta*r12));

        }
    }

    return exp(wf);
}


double correlation::getRatioJastrow(int i, const mat &rOld, const mat &rNew, double beta) {

    double exponent=0;
    double a=0;
    double rijOld = 0;
    double rijNew = 0;

    int test1 = i-nParticles/2;
    for(int j=0;j<nParticles;j++) {
        if(j != i) {
            int test2 = j-nParticles/2;
            if(test1*test2>0) a=0.25;
            else a=0.5;
            rijNew = 0; rijOld = 0;
            for(int d=0;d<nDimensions;d++) {
                rijNew += pow((rNew(i,d)-rNew(j,d)),2);
                rijOld += pow((rOld(i,d)-rOld(j,d)),2);
            }
            rijNew = sqrt(rijNew);
            rijOld = sqrt(rijOld);
            exponent += (a*rijNew)/(1+beta*rijNew) - (a*rijOld)/(1+beta*rijOld);

        }
    }


    return exp(exponent);

}


vec correlation::gradientWaveFunction(const mat &r, int i, double beta) {

    vec grad = zeros<vec>(nDimensions);
    double a = 0.25;
    int test1 = i-nParticles/2;
    double rij = 0;

    for(int j=0; j<nParticles; j++) {
        if(j!=i) {
            int test2 = j-nParticles/2;
            if(test1*test2>0) a=0.25;
            else a=0.5;
            rij = 0;
            for (int d=0; d<nDimensions; d++) { rij += pow((r(i,d)-r(j,d)),2); } //Get r for particle
            rij = sqrt(rij);
            grad(0) += ((r(i,0)-r(j,0))/rij)*(a/pow((1+beta*rij),2));
            grad(1) += ((r(i,1)-r(j,1))/rij)*(a/pow((1+beta*rij),2));
            grad(2) += ((r(i,2)-r(j,2))/rij)*(a/pow((1+beta*rij),2));
        }
    }

    return grad;

}

vec correlation::gradientWaveFunctionNum(const mat &r, int i, double beta) {

    vec grad = zeros(nDimensions);
    mat rPlus = zeros<mat>(nDimensions);
    mat rMinus = zeros<mat>(nDimensions);
    double h = 0.001;

    double wf = jastrowNum(r,beta);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    //First derivative

    for(int p=0; p<nParticles; p++) {
        if(p != i) {
            for(int j = 0; j < nDimensions; j++) {
                rPlus(i,j) = r(i,j)+h;
                rMinus(i,j) = r(i,j)-h;
                waveFunctionMinus = jastrowNum(rMinus,beta);
                waveFunctionPlus = jastrowNum(rPlus,beta);
                grad(j) = (waveFunctionPlus - waveFunctionMinus)/(2*wf*h);
                rPlus(i,j) = r(i,j);
                rMinus(i,j) = r(i,j);
            }
        }
    }
    //cout <<"Numerical: "<< grad << endl;

    return grad;
}



double correlation::laPlaceWaveFunction(const mat &r, double beta) {

    double rik = 0;
    double rjk = 0;
    double laplace = 0;
    double a1 = 0.25;
    double a2 = 0.25;
    int test1;
    int test2;
    int test3;
    vec rik_vec = zeros<vec>(nDimensions,1);
    vec rjk_vec = zeros<vec>(nDimensions,1);

    for(int k=0; k<nParticles; k++) {
        test1 = k-nParticles/2;
        for(int i = 0; i < nParticles; i++) {
            test2 = i-nParticles/2;
            if(test1*test2>0) a1=0.25;
            else a1=0.5;
            for(int j=0; j<nParticles; j++) {
                if(i != k && j != k) {
                    test3 = j-nParticles/2;
                    if(test1*test3>0) a2=0.25;
                    else a2=0.5;
                    rik = 0; rjk = 0;
                    for(int d=0;d<nDimensions;d++) {
                        rik += pow((r(i,d)-r(k,d)),2);
                        rjk += pow((r(j,d)-r(k,d)),2);
                        rik_vec(d) = r(k,d) - r(i,d);
                        rjk_vec(d) = r(k,d) - r(j,d);
                    }
                    rik = sqrt(rik);
                    rjk = sqrt(rjk);
                    laplace += (dot(rik_vec,rjk_vec)/(rik*rjk))*(a1/pow((1+beta*rik),2)*(a2/pow((1+beta*rjk),2)));
                }
            }
        }
    }

    for(int k=0; k<nParticles; k++) {
        test1 = k-nParticles/2;
        for(int j=0; j<nParticles; j++) {
            if(j != k) {
                int test3 = j-nParticles/2;
                if(test1*test3>0) a2=0.25;
                else a2=0.5;
                rjk = 0;
                for(int d=0;d<nDimensions;d++) {
                    rjk += pow((r(j,d)-r(k,d)),2);
                    rjk_vec(d) = r(k,d) - r(j,d);
                }
                rjk = sqrt(rjk);
                laplace += (2*a2)/(rjk*pow((1+beta*rjk),2)) - (2*a2*beta)/pow((1+beta*rjk),3);
            }
        }
    }

    return laplace;

}

double correlation::laPlaceWaveFunctionNum(const mat &r, double beta) {

    double h2 = 1000000;
    double h = 0.001;

    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    double waveFunctionCurrent = jastrowNum(r,beta);

    //Second derivative (del^2):

    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = jastrowNum(rMinus,beta);
            waveFunctionPlus = jastrowNum(rPlus,beta);

           kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
    kineticEnergy =  h2 * kineticEnergy / waveFunctionCurrent;

    //cout <<"Numerical: "<<kineticEnergy<<endl;

    return kineticEnergy;


}




