#include "vmcsolver.h"
#include "lib.h"
#include "WaveFunction.h"
#include "hamiltonian.h"
#include "slaterdeterminant.h"
#include "correlation.h"
#include "minimise.h"


#include <armadillo>
#include <fstream>
#include <iostream>
#include <mpi.h>


using namespace arma;
using namespace std;

VMCSolver::VMCSolver():
    nDimensions(3),
    h(0.001),
    h2(1000000),
    timestep(0.05),
    D(0.5),
    stepLength(1.0),
    nCycles(100000),
    charge(4),
    nParticles(4),
    alpha(4),
    beta(0.5),
    minimise_var(true),
    min_steps(10000)     //Number of steps for minimiser


{
}

void VMCSolver::runMonteCarloIntegration(int argc, char *argv[])
{
    bool printToFile = false;

    char file_energies[] = "../../../output/energy.txt";
    char file_energySquareds[] = "../../../output/squareds.txt";
    char file_alpha[] = "../../../output/alpha_beta.txt";
    char file_sigma[] = "../../../output/sigma.txt";

    slaterDeterminant *slater = new slaterDeterminant(nParticles, nDimensions);
    Hamiltonian *hamiltonian = new Hamiltonian(nParticles, nDimensions, h, h2, charge);
    correlation *corr = new correlation(nParticles, nDimensions);

    double energies = 0;
    double energySquareds = 0;
    long idum = -1;
    int id, np;


    //mpd --ncpus=4 &
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    double myTime,mintime, maxtime,avgtime;
    myTime = MPI_Wtime();

    int mpi_steps = nCycles/np;
    idum = idum-id;

    double* allEnergies = new double[mpi_steps+1];
    double pEnergies = 0;
    double pEnergySquareds = 0;
    vec energySums = zeros<vec>(2);

    cout << "ID: " << id << endl;


    if(minimise_var) {

        double gtol = 1e-3;
        int iter;
        double fret;
        vec p = zeros<vec>(2,1);
        p(0) = alpha;
        p(1) = beta;
        int n = 2;

        //dfpmin(idum, p, n, gtol, min_steps, &iter, &fret, slater,corr, hamiltonian);
        //cout <<"E, iter: "<<fret<<" "<<iter<<endl;
        //cout <<p<<endl;

        vec ans = steepest_descent(idum, p, n, gtol, min_steps, &iter, &fret, slater,corr, hamiltonian);
        cout <<ans<<endl;
        double alpha_new = ans(0);
        double beta_new = ans(1);

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Allreduce(&alpha_new, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&beta_new, &beta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        alpha = alpha/np;
        beta = beta/np;
        cout << "Final alpha, beta: "<< alpha<<" "<<beta<<endl;

//        alpha = p(0);
//        beta = p(1);
    }

    energySums = MCImportance(idum, alpha, beta, mpi_steps, slater, hamiltonian, corr, allEnergies);
    //energySums = MCSampling(idum, alpha, beta, mpi_steps, slater, hamiltonian, corr, allEnergies);

    if(printToFile) {
        ostringstream ost;
        //ost << "/mn/korona/rp-s1/alborg/4411/helium/examples/vmc-simple/DATA/data" << id << ".mat" ;
        ost << "../vmc-simple/DATA/data" << id << ".mat" ;

        ofstream blockofile;
        blockofile.open( ost.str( ).c_str( ),ios::out | ios::binary );
        if (blockofile.is_open())
        {
            blockofile.write((char*)(allEnergies+1) , mpi_steps*sizeof(double)) ;
            blockofile.close();
        }
        else cout << "Unable to open data file for process " << id << endl;
    }

    pEnergies = energySums(0)/(nCycles * nParticles);
    pEnergySquareds = energySums(1)/(nCycles * nParticles);

    cout << "--------------------------" << endl;



    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allreduce(&pEnergies, &energies, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&pEnergySquareds, &energySquareds, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    myTime = MPI_Wtime() - myTime;
    MPI_Reduce(&myTime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&myTime, &mintime, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
    MPI_Reduce(&myTime, &avgtime, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
    MPI_Finalize();

    if (id == 0) {
        cout << "Energies: " << energies << endl; //*2*13.6
        cout << "Energy squareds: " << energySquareds << endl; //*2*13.6*2*13.6
       // printFile(*file_energies, *file_energySquareds, *file_alpha, *file_sigma, energies, energySquareds, alphas, betas);
        avgtime /= np;
        cout << "Min time: " << mintime << ", max time: " << maxtime << ", avg time: " << avgtime << endl;
    }

    delete[] allEnergies;
}

vec VMCSolver::MCImportance(long idum, double alpha, double beta, int mpi_steps,
                             slaterDeterminant *slater, Hamiltonian *hamiltonian, correlation *corr,
                             double *allEnergies) {



    vec qForceOld = zeros<vec>(nDimensions,1);
    vec qForceNew = zeros<vec>(nDimensions,1);
    mat rOld = zeros<mat>(nParticles, nDimensions);
    mat rNew = zeros<mat>(nParticles, nDimensions);
    double accepted_steps = 0;
    double count_total = 0;
    double deltaE = 0;
    vec deltaPsi = zeros<vec>(2);
    vec deltaPsiE = zeros<vec>(2);
    double cycleE = 0;
    double ratio = 1;
    double ratioCorr = 1;
    vec energySums = zeros<vec>(6);

    // initial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = gaussianDeviate(&idum)*sqrt(timestep);
        }
    }


    //Build the full Slater matrix (done only once):
    slater->buildDeterminant(rOld, alpha, beta);

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < mpi_steps; cycle++) {


        for(int i = 0; i < nParticles; i++) { //Particle

            ratio = 1.0;
            qForceOld = slater->gradientWaveFunction(rOld, i, ratio, alpha, beta);

            // New position to test
            for(int d = 0; d < nDimensions; d++) {
                rNew(i,d) = rOld(i,d) + gaussianDeviate(&idum)*sqrt(timestep) + qForceOld(d)*timestep*D;
            }

            //Move only one particle (i).
            for (int g=0; g<nParticles; g++) {
                if(g != i) {
                    for(int d=0; d<nDimensions; d++) {
                        rNew(g,d) = rOld(g,d);
                    }
                }
            }

            //Get the ratio of the new to the old determinant (wavefunction).
            ratio = slater->getRatioDeterminant(i, rNew, alpha, beta);
            ratioCorr = corr->getRatioJastrow(i, rOld, rNew, beta);

            qForceNew = slater->gradientWaveFunction(rNew, i, ratio, alpha, beta);

            //Greens function
            double greensFunction = 0;
            for(int d=0; d<nDimensions; d++) {
                greensFunction += 0.5*(qForceOld(d) + qForceNew(d)) * (0.5*D*timestep*(qForceOld(d) - qForceNew(d)) - rNew(i,d) + rOld(i,d));
            }
            greensFunction = exp(greensFunction);

            ++count_total;


            // Check for step acceptance (if yes, update position and determinant, if no, reset position)
           if(ran2(&idum) <= greensFunction * ratio*ratio * ratioCorr*ratioCorr) {
                ++accepted_steps;
               slater->updateDeterminant(rNew, rOld, i, alpha, beta, ratio);
               for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);   
                }
            }
            else {
               for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j);
                }
            }


            // update energies
            deltaE = hamiltonian->localEnergy(rNew, alpha, beta, slater,corr);
            energySums(0) += deltaE;
            energySums(1) += deltaE*deltaE;
            allEnergies[cycle] += deltaE;
            cycleE += deltaE;
            if(minimise_var) {
                deltaPsi = hamiltonian->dPsi(rNew,alpha,beta,slater,corr);
                deltaPsiE(0) = deltaE*deltaPsi(0);
                deltaPsiE(1) = deltaE*deltaPsi(1);
                energySums(2) += deltaPsi(0);
                energySums(3) += deltaPsi(1);
                energySums(4) += deltaPsiE(0);
                energySums(5) += deltaPsiE(1);
            }


        } //End particle loop

        allEnergies[cycle] += cycleE;
        cycleE = 0;

    } //End Monte Carlo loop

    cout << "accepted steps: " << 100*accepted_steps/count_total << "%" << endl;

    return energySums;

}


vec VMCSolver::MCSampling(long idum, double alpha, double beta, int mpi_steps,
                           slaterDeterminant *slater, Hamiltonian *hamiltonian,
                           correlation *corr, double *allEnergies) {

    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);
    int accepted_steps = 0;
    int count_total = 0;
    double deltaE = 0;
    double cycleE = 0;
    vec energySums = zeros<vec>(2);



    // initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
        }
    }
    rNew = rOld;

     slater->buildDeterminant(rOld, alpha, beta);

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < mpi_steps; cycle++) {


        // New position to test
        for(int i = 0; i < nParticles; i++) { //Particles

            for(int d = 0; d < nDimensions; d++) {
                rNew(i,d) = rOld(i,d) + stepLength*(ran2(&idum) - 0.5);
            }

            ++count_total;

            double ratioSlater = slater->getRatioDeterminant(i, rNew, alpha, beta);
            double ratioCorr = corr->getRatioJastrow(i, rOld, rNew, beta);

            // Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= ratioSlater*ratioSlater*ratioCorr*ratioCorr) {
                ++accepted_steps;
                slater->updateDeterminant(rNew, rOld, i, alpha, beta, ratioSlater);
                for(int d = 0; d < nDimensions; d++) {
                    rOld(i,d) = rNew(i,d);
                }
            }
            else {
                for(int d = 0; d < nDimensions; d++) {
                    rNew(i,d) = rOld(i,d);
                }
            }


            // update energies
            deltaE = hamiltonian->localEnergy(rNew, alpha, beta, slater,corr);
            //deltaE = hamiltonian->analyticEnergyH(rNew, alpha, beta);
            energySums(0) += deltaE;
            energySums(1) += deltaE*deltaE;
            cycleE += deltaE;

        } //Particles

        allEnergies[cycle] += cycleE;
        cycleE = 0;

    } //Monte Carlo cycles

    cout << "accepted steps: " << 100*accepted_steps/count_total << "%" << endl;

    return energySums;
}





double VMCSolver::gaussianDeviate(long *idum) {

    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;

    if ( idum < 0) iset =0;
    if (iset == 0) {
        do {
            v1 = 2.*ran2(idum) -1.0;
            v2 = 2.*ran2(idum) -1.0;
            rsq = v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0.);
        fac = sqrt(-2.*log(rsq)/rsq);
        gset = v1*fac;
        iset = 1;
        return v2*fac;
    } else {
        iset =0;
        return gset;
    }


}

vec VMCSolver::steepest_descent(long idum, vec &p, int n, double gtol, int min_steps, int *iter, double *fret,
                       slaterDeterminant *slater, correlation *corr, Hamiltonian *hamiltonian)
{
    vec dPsi = zeros<vec>(2,1);
    vec dPsi_Elocal = zeros<vec>(2,1);
    double* allEnergies = new double[min_steps+1];
    double alpha = p(0);
    double beta = p(1);
    double alpha_new = alpha;
    double beta_new = beta;
    vec dE = zeros<vec>(n);
    vec dEold = zeros<vec>(n);
    int maxIter = 20;
    vec answers = zeros<vec>(n+2);
    double E = 0;
    double Enew = 0;
    double alpha_step = 0.5;
    double beta_step = 0.5;
    double step = 0.5;
    int i = 0;
    double test;
    double step_reduce = 1.5;

    cout<<"Steepest"<<endl;

    vec Es = MCImportance(idum, alpha,beta,min_steps, slater, hamiltonian, corr, allEnergies);
    E = Es(0)/(min_steps * nParticles);
    dPsi(0) = Es(2)/(min_steps * nParticles);
    dPsi(1) = Es(3)/(min_steps * nParticles);
    dPsi_Elocal(0) = Es(4)/(min_steps * nParticles);
    dPsi_Elocal(1) = Es(5)/(min_steps * nParticles);
    gradE(dPsi, E, dPsi_Elocal,dE);

    cout <<"E: "<<E<<endl;

    while(i<maxIter) {

        //Alpha:

        alpha_new = alpha - step*dE(0);
        if(alpha_new < 0) {
//            while(alpha_new < 0) {
//                step = step/2;
//                alpha_new = alpha - step*dE(0);
//            }
            alpha_new = alpha;
        }

        Es = MCImportance(idum, alpha_new,beta,min_steps, slater, hamiltonian, corr, allEnergies);
        Enew = Es(0)/(min_steps * nParticles);
        dEold(0) = dE(0);
        dPsi(0) = Es(2)/(min_steps * nParticles);
        dPsi(1) = Es(3)/(min_steps * nParticles);
        dPsi_Elocal(0) = Es(4)/(min_steps * nParticles);
        dPsi_Elocal(1) = Es(5)/(min_steps * nParticles);
        gradE(dPsi, E, dPsi_Elocal,dE);
        if(dE(0)*dEold(0) < 0) step = step/step_reduce;

        cout <<"Alpha new: "<< alpha_new << " " << beta<<" "<<endl;
        cout<<"Enew: "<<Enew<<" "<<step<<endl;

        test = abs(Enew-E);
        if(test < 1e-5) break;

        if(Enew < E) {
            E = Enew;
            //if(dE(0)*dEold(0) > 0) step = step/step_reduce;
        }

        alpha = alpha_new;

        //Beta:
        beta_new = beta - step*dE(1);
        if(beta_new < 0) {
            beta_new = beta;

        }


        Es = MCImportance(idum, alpha,beta_new,min_steps, slater, hamiltonian, corr, allEnergies);
        Enew = Es(0)/(min_steps * nParticles);
        dEold(1) = dE(1);
        dPsi(0) = Es(2)/(min_steps * nParticles);
        dPsi(1) = Es(3)/(min_steps * nParticles);
        dPsi_Elocal(0) = Es(4)/(min_steps * nParticles);
        dPsi_Elocal(1) = Es(5)/(min_steps * nParticles);
        gradE(dPsi, E, dPsi_Elocal,dE);
        if(dE(1)*dEold(1) < 0) step = step/step_reduce;

        cout <<"Beta new: "<< alpha<<" "<<beta_new<<" "<<endl;
        cout<<"Enew: "<<Enew<<" "<<step<<endl;

        test = abs(Enew-E);
        if(test < 1e-5) break;

        if(Enew < E) {
            E = Enew;
            //if(dE(1)*dEold(1) > 0) step = step/step_reduce;

        }

        beta = beta_new;

        i++;
    }

    answers(0) = alpha_new;
    answers(1) = beta_new;
    answers(2) = Enew;
    answers(3) = i;

    return answers;
}


double VMCSolver::gradE(vec dPsi, double Elocal, vec dPsi_Elocal, vec &g) {

    g(0) = 2*(dPsi(0)*Elocal - dPsi_Elocal(0));
    g(1) = 2*(dPsi(1)*Elocal - dPsi_Elocal(1));
}



static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))



#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0


void VMCSolver::dfpmin(long idum, vec &p, int n, double gtol, int min_steps, int *iter, double *fret,
                       slaterDeterminant *slater, correlation *corr, Hamiltonian *hamiltonian)
{
    int check,i,its,j;
    double den,fac,fad,fp,fae,stpmax,sum=0.0,sumdg,sumxi,temp,test;
    vec dg = zeros<vec>(n);
    vec g = zeros<vec>(n);
    vec hdg = zeros<vec>(n);
    vec pnew = zeros<vec>(n);
    vec xi = zeros<vec>(n);
    mat hessian = zeros<vec>(n,n);
    vec dPsi = zeros<vec>(2,1);
    vec dPsi_Elocal = zeros<vec>(2,1);
    double* allEnergies = new double[min_steps+1];
    double alpha_start = p(0);
    double beta_start = p(1);

    vec Es = MCImportance(idum, p(0),p(1),min_steps, slater, hamiltonian, corr, allEnergies);
    fp = Es(0)/(min_steps * nParticles);
    dPsi(0) = Es(2)/(min_steps * nParticles);
    dPsi(1) = Es(3)/(min_steps * nParticles);
    dPsi_Elocal(0) = Es(4)/(min_steps * nParticles);
    dPsi_Elocal(1) = Es(5)/(min_steps * nParticles);
    gradE(dPsi, fp, dPsi_Elocal,g);

    for (i = 0;i < n;i++) {
        for (j = 0; j< n;j++) hessian(i,j)=0.0;
        hessian(i,i)=1.0;
        xi(i) = -g(i);
        sum += p(i)*p(i);
    }
    stpmax=STPMX*FMAX(sqrt(sum),(double)n);
    for (its=1;its<=ITMAX;its++) {
        *iter=its;
        lnsrch(idum, alpha_start, beta_start,n,p,fp,g,xi,pnew,fret,stpmax,&check,min_steps, slater, hamiltonian, corr, allEnergies);
        fp = *fret;
        for (i = 0; i< n;i++) {
            xi(i)=pnew(i)-p(i);
            p(i)=pnew(i);
        }
        if(p(0) < 0) p(0) = alpha_start;
        if(p(1) < 0) p(1) = beta_start;
        test=0.0;
        for (i = 0;i< n;i++) {
            temp=fabs(xi(i))/FMAX(fabs(p(i)),1.0);
            if (temp > test) test=temp;
        }
        cout<<"alpha, beta 1: "<< p <<endl;
        if (test < TOLX) {
            return;
        }
        for (i=0;i<n;i++) dg(i)=g(i);
        gradE(dPsi, fp, dPsi_Elocal,g);
        test=0.0;
        den=FMAX(*fret,1.0);
        for (i=0;i<n;i++) {
            temp=fabs(g(i))*FMAX(fabs(p(i)),1.0)/den;
            if (temp > test) test=temp;
        }
        cout<<"alpha, beta 2: "<< p <<endl;
        if (test < gtol) {
            return;
        }
        for (i=0;i<n;i++) dg(i)=g(i)-dg(i);
        for (i=0;i<n;i++) {
            hdg(i)=0.0;
            for (j=0;j<n;j++) hdg(i) += hessian(i,j)*dg(j);
        }
        fac=fae=sumdg=sumxi=0.0;
        for (i=0;i<n;i++) {
            fac += dg(i)*xi(i);
            fae += dg(i)*hdg(i);
            sumdg += SQR(dg(i));
            sumxi += SQR(xi(i));
        }
        if (fac*fac > EPS*sumdg*sumxi) {
            fac=1.0/fac;
            fad=1.0/fae;
            for (i=0;i<n;i++) dg(i)=fac*xi(i)-fad*hdg(i);
            for (i=0;i<n;i++) {
                for (j=0;j<n;j++) {
                    hessian(i,j) += fac*xi(i)*xi(j)
                            -fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j);
                }
            }
        }
        for (i=0;i<n;i++) {
            xi(i)=0.0;
            for (j=0;j<n;j++) xi(i) -= hessian(i,j)*g(j);
        }
    }
    cout << "too many iterations in dfpmin" << endl;
}
#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX

#define ALF 1.0e-4
#define TOLX 1.0e-7

void VMCSolver::lnsrch(long idum, double alpha_start, double beta_start, int n, vec &xold,double fold,
                       vec &g,vec &p,vec &x,double *f, double stpmax,
                       int *check,int min_steps, slaterDeterminant *slater, Hamiltonian *hamiltonian,
                       correlation *corr, double *allEnergies)
{

    int i;
    double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
            test,tmplam;

    *check=0;
    for (sum=0.0,i=0;i<n;i++) sum += p(i)*p(i);
    sum=sqrt(sum);
    if (sum > stpmax)
        for (i=0;i<n;i++) p(i) *= stpmax/sum;
    for (slope=0.0,i=0;i<n;i++)
        slope += g(i)*p(i);
    test=0.0;
    for (i=0;i<n;i++) {
        temp=fabs(p(i))/FMAX(fabs(xold(i)),1.0);
        if (temp > test) test=temp;
    }
    alamin=TOLX/test;
    alam=1.0;
    for (;;) {
        for (i=0;i<n;i++) {
            x(i)=xold(i)+alam*p(i);
        }
        if(x(0) < 0) x(0) = alpha_start;
        if(x(1) < 0) x(1) = beta_start;
        cout << x << endl;
        vec Es = MCImportance(idum, x(0),x(1),min_steps, slater, hamiltonian, corr, allEnergies);
        *f = Es(0)/(min_steps * nParticles);
        if (alam < alamin) {
            for (i=0;i<n;i++) x(i)=xold(i);
            *check=1;
            return;
        } else if (*f <= fold+ALF*alam*slope) return;
        else {
            if (alam == 1.0)
                tmplam = -slope/(2.0*(*f-fold-slope));
            else {
                rhs1 = *f-fold-alam*slope;
                rhs2=f2-fold2-alam2*slope;
                a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
                b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
                if (a == 0.0) tmplam = -slope/(2.0*b);
                else {
                    disc=b*b-3.0*a*slope;
                    if (disc<0.0) cout << "Roundoff problem in lnsrch." << endl;
                    else tmplam=(-b+sqrt(disc))/(3.0*a);
                }
                if (tmplam>0.5*alam)
                    tmplam=0.5*alam;
            }
        }
        alam2=alam;
        f2 = *f;
        fold2=fold;
        alam=FMAX(tmplam,0.1*alam);
    }
}
#undef ALF
#undef TOLX


void VMCSolver::printFile(const char &file_energies, const char &file_energySquareds, const char &file_sigma, const char &file_alpha, const mat &energies, const mat &energiesSquared, const vec alphas, const vec betas)
{

    ofstream myfile(&file_energies);
    if (myfile.is_open())
    {
        for (unsigned int f=0; f<energies.n_rows; f++)
        {
            for (unsigned int l=0; l<energies.n_cols; l++) {
                myfile << energies(f,l) << " ";
            }
            myfile << endl;
        }

        myfile.close();
    }
    else cout << "Unable to open file" << endl;


    ofstream myfile2 (&file_alpha);
    if (myfile2.is_open())
    {
        myfile2 << alphas << endl;
        myfile2 << betas << endl;

        myfile2.close();
    }
    else cout << "Unable to open file" << endl;


    ofstream myfile3(&file_energySquareds);
    if (myfile3.is_open())
    {
        for (unsigned int f=0; f<energiesSquared.n_rows; f++)
        {
            for (unsigned int l=0; l<energiesSquared.n_cols; l++) {
                myfile3 << energiesSquared(f,l) << " ";
            }
            myfile3 << endl;
        }

        myfile3.close();
    }
    else cout << "Unable to open file" << endl;


    ofstream myfile4(&file_sigma);
    if (myfile4.is_open())
    {
        for (unsigned int f=0; f<energiesSquared.n_rows; f++)
        {
            for (unsigned int l=0; l<energiesSquared.n_cols; l++) {
                myfile4 << sqrt(energiesSquared(f,l) - energies(f,l)*energies(f,l))<< " ";
            }
            myfile4 << endl;
        }

        myfile4.close();
    }
    else cout << "Unable to open file" << endl;

}

