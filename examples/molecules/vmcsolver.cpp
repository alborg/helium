#include "vmcsolver.h"
#include "lib.h"
#include "WaveFunction.h"
#include "hamiltonian.h"


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
    nCycles(10000),
    alpha(1.3),
    beta(0.3),
    R(1.34),
    nProtons(2),
    nElectrons(1),
    nParticles(nProtons*nElectrons),
    charge(1),
    minimise_var(false),
    min_steps(1000),     //Number of steps for minimiser
    printToFile(false)

{
}

void VMCSolver::runMonteCarloIntegration(int argc, char *argv[])
{

    WaveFunction *function = new WaveFunction(nDimensions,nProtons,nElectrons);
    Hamiltonian *hamiltonian = new Hamiltonian(nProtons, nElectrons,nDimensions,h,h2,charge);

    double energies = 0;
    double energySquareds = 0;
    long idum = -1;

    int id, np;

    rProtons = zeros<mat>(nProtons, nDimensions);
    rProtons(0,2) = -R/2;
    rProtons(1,2) = R/2;

//    if(minimise_var) {

//        double gtol = 1e-5;

//        int iter;
//        double fret;
//        vec p = zeros<vec>(2,1);
//        p(0) = alpha;
//        p(1) = beta;
//        int n = 2;

//        dfpmin(p, n, gtol, min_steps, &iter, &fret, hamiltonian);

//        cout <<"E, iter: "<<fret<<" "<<iter<<endl;
//        cout <<p<<endl;
//    }


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

    energySums = MCImportance(idum, alpha, beta, mpi_steps, hamiltonian, function, allEnergies);

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
        avgtime /= np;
        cout << "Min time: " << mintime << ", max time: " << maxtime << ", avg time: " << avgtime << endl;
    }

    delete[] allEnergies;
}

vec VMCSolver::MCImportance(long idum, double alpha, double beta, int mpi_steps, Hamiltonian *hamiltonian,
                            WaveFunction *function, double *allEnergies) {


    mat rOld = zeros<mat>(nParticles, nDimensions);
    mat rNew = zeros<mat>(nParticles, nDimensions);
    double accepted_steps = 0;
    double count_total = 0;
    double deltaE = 0;
//    vec deltaPsi = zeros<vec>(2);
//    vec deltaPsiE = zeros<vec>(2);
    double cycleE = 0;
    vec energySums = zeros<vec>(6);

    mat qForceOld = zeros(nParticles,nDimensions);
    mat qForceNew = zeros(nParticles,nDimensions);
    double waveFunctionOld = 0;
    double waveFunctionNew = 0;


    // initial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = gaussianDeviate(&idum)*sqrt(timestep);
        }
    }

    waveFunctionOld = function->waveFunction(rOld, rProtons, alpha, beta);
    qForceOld = quantumForce(rOld, rProtons, alpha, beta, waveFunctionOld,function);

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < mpi_steps; cycle++) {


        for(int i = 0; i < nParticles; i++) { //Particle

            // New position to test
            for(int d = 0; d < nDimensions; d++) {
                rNew(i,d) = rOld(i,d) + gaussianDeviate(&idum)*sqrt(timestep) + qForceOld(i,d)*timestep*D;
            }

            //Move only one particle (i).
            for (int g=0; g<nParticles; g++) {
                if(g != i) {
                    for(int d=0; d<nDimensions; d++) {
                        rNew(g,d) = rOld(g,d);
                    }
                }
            }

            // Recalculate the wave function
            waveFunctionNew = function->waveFunction(rNew, rProtons, alpha, beta);
            qForceNew = quantumForce(rNew, rProtons, alpha, beta, waveFunctionNew,function);


            //Greens function
            double greensFunction = 0;
            for(int d=0; d<nDimensions; d++) {
                greensFunction += 0.5*(qForceOld(i,d) + qForceNew(i,d)) * (0.5*D*timestep*(qForceOld(i,d) - qForceNew(i,d)) - rNew(i,d) + rOld(i,d));
            }
            greensFunction = exp(greensFunction);

            ++count_total;


            // Check for step acceptance (if yes, update position and determinant, if no, reset position)
           if(ran2(&idum) <= greensFunction * (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                ++accepted_steps;
               for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                    qForceOld(i,j) = qForceNew(i,j);
                }
               waveFunctionOld = waveFunctionNew;
            }
            else {
               for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j);
                    qForceNew(i,j) = qForceOld(i,j);
                }
            }


            // update energies
            deltaE = hamiltonian->localEnergy(R, rNew, rProtons, alpha, beta, function);
            energySums(0) += deltaE;
            energySums(1) += deltaE*deltaE;
            allEnergies[cycle] += deltaE;
            cycleE += deltaE;
            if(minimise_var) {
                //deltaPsi = hamiltonian->dPsi(rNew,alpha,beta,slater,corr);
//                deltaPsiE(0) = deltaE*deltaPsi(0);
//                deltaPsiE(1) = deltaE*deltaPsi(1);
//                energySums(2) += deltaPsi(0);
//                energySums(3) += deltaPsi(1);
//                energySums(4) += deltaPsiE(0);
//                energySums(5) += deltaPsiE(1);
            }


        } //End particle loop

        allEnergies[cycle] += cycleE;
        cycleE = 0;

    } //End Monte Carlo loop

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


mat VMCSolver::quantumForce(const mat &r, const mat &rProtons, double alpha_, double beta_, double wf, WaveFunction *function) {

    mat qforce = zeros(nParticles, nDimensions);
    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    //First derivative

    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) = r(i,j)+h;
            rMinus(i,j) = r(i,j)-h;
            waveFunctionMinus = function->waveFunction(rMinus, rProtons, alpha_, beta_);
            waveFunctionPlus = function->waveFunction(rPlus, rProtons, alpha_, beta_);
            qforce(i,j) = (waveFunctionPlus - waveFunctionMinus)/(wf*h);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }

    return qforce;
}



//double VMCSolver::gradE(vec dPsi, double Elocal, vec dPsi_Elocal, vec &g) {

//    g(0) = 2*(dPsi(0)*Elocal - dPsi_Elocal(0));
//    g(1) = 2*(dPsi(1)*Elocal - dPsi_Elocal(1));
//}



//static double sqrarg;
//#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


//static double maxarg1,maxarg2;
//#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
//        (maxarg1) : (maxarg2))



//#define ITMAX 200
//#define EPS 3.0e-8
//#define TOLX (4*EPS)
//#define STPMX 100.0


//void VMCSolver::dfpmin(vec &p, int n, double gtol, int min_steps, int *iter, double *fret, Hamiltonian *hamiltonian)
//{
//    int check,i,its,j;
//    double den,fac,fad,fp,fae,stpmax,sum=0.0,sumdg,sumxi,temp,test;
//    vec dg = zeros<vec>(n);
//    vec g = zeros<vec>(n);
//    vec hdg = zeros<vec>(n);
//    vec pnew = zeros<vec>(n);
//    vec xi = zeros<vec>(n);
//    mat hessian = zeros<vec>(n,n);
//    vec dPsi = zeros<vec>(2,1);
//    vec dPsi_Elocal = zeros<vec>(2,1);
//    double* allEnergies = new double[min_steps+1];
//    double alpha_start = p(0);
//    double beta_start = p(1);
//    long idum = -1;

//    vec Es = MCImportance(idum, p(0),p(1),min_steps, slater, hamiltonian, corr, allEnergies);
//    fp = Es(0)/(min_steps * nParticles);
//    dPsi(0) = Es(2)/(min_steps * nParticles);
//    dPsi(1) = Es(3)/(min_steps * nParticles);
//    dPsi_Elocal(0) = Es(4)/(min_steps * nParticles);
//    dPsi_Elocal(1) = Es(5)/(min_steps * nParticles);
//    gradE(dPsi, fp, dPsi_Elocal,g);

//    for (i = 0;i < n;i++) {
//        for (j = 0; j< n;j++) hessian(i,j)=0.0;
//        hessian(i,i)=1.0;
//        xi(i) = -g(i);
//        sum += p(i)*p(i);
//    }
//    stpmax=STPMX*FMAX(sqrt(sum),(double)n);
//    for (its=1;its<=ITMAX;its++) {
//        *iter=its;
//        idum -= 1;
//        lnsrch(idum, alpha_start, beta_start,n,p,fp,g,xi,pnew,fret,stpmax,&check,min_steps, slater, hamiltonian, corr, allEnergies);
//        fp = *fret;
//        for (i = 0; i< n;i++) {
//            xi(i)=pnew(i)-p(i);
//            p(i)=pnew(i);
//        }
//        if(p(0) < 0) p(0) = alpha_start;
//        if(p(1) < 0) p(1) = beta_start;
//        test=0.0;
//        for (i = 0;i< n;i++) {
//            temp=fabs(xi(i))/FMAX(fabs(p(i)),1.0);
//            if (temp > test) test=temp;
//        }
//        cout<<"alpha, beta 1: "<< p <<endl;
//        if (test < TOLX) {
//            return;
//        }
//        for (i=0;i<n;i++) dg(i)=g(i);
//        gradE(dPsi, fp, dPsi_Elocal,g);
//        test=0.0;
//        den=FMAX(*fret,1.0);
//        for (i=0;i<n;i++) {
//            temp=fabs(g(i))*FMAX(fabs(p(i)),1.0)/den;
//            if (temp > test) test=temp;
//        }
//        cout<<"alpha, beta 2: "<< p <<endl;
//        if (test < gtol) {
//            return;
//        }
//        for (i=0;i<n;i++) dg(i)=g(i)-dg(i);
//        for (i=0;i<n;i++) {
//            hdg(i)=0.0;
//            for (j=0;j<n;j++) hdg(i) += hessian(i,j)*dg(j);
//        }
//        fac=fae=sumdg=sumxi=0.0;
//        for (i=0;i<n;i++) {
//            fac += dg(i)*xi(i);
//            fae += dg(i)*hdg(i);
//            sumdg += SQR(dg(i));
//            sumxi += SQR(xi(i));
//        }
//        if (fac*fac > EPS*sumdg*sumxi) {
//            fac=1.0/fac;
//            fad=1.0/fae;
//            for (i=0;i<n;i++) dg(i)=fac*xi(i)-fad*hdg(i);
//            for (i=0;i<n;i++) {
//                for (j=0;j<n;j++) {
//                    hessian(i,j) += fac*xi(i)*xi(j)
//                            -fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j);
//                }
//            }
//        }
//        for (i=0;i<n;i++) {
//            xi(i)=0.0;
//            for (j=0;j<n;j++) xi(i) -= hessian(i,j)*g(j);
//        }
//    }
//    cout << "too many iterations in dfpmin" << endl;
//}
//#undef ITMAX
//#undef EPS
//#undef TOLX
//#undef STPMX

//#define ALF 1.0e-4
//#define TOLX 1.0e-7

//void VMCSolver::lnsrch(long idum, double alpha_start, double beta_start, int n, vec &xold, double fold,
//                       vec &g, vec &p, vec &x, double *f, double stpmax,
//                       int *check, int min_steps, Hamiltonian *hamiltonian, double *allEnergies)
//{

//    int i;
//    double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
//            test,tmplam;

//    *check=0;
//    for (sum=0.0,i=0;i<n;i++) sum += p(i)*p(i);
//    sum=sqrt(sum);
//    if (sum > stpmax)
//        for (i=0;i<n;i++) p(i) *= stpmax/sum;
//    for (slope=0.0,i=0;i<n;i++)
//        slope += g(i)*p(i);
//    test=0.0;
//    for (i=0;i<n;i++) {
//        temp=fabs(p(i))/FMAX(fabs(xold(i)),1.0);
//        if (temp > test) test=temp;
//    }
//    alamin=TOLX/test;
//    alam=1.0;
//    for (;;) {
//        for (i=0;i<n;i++) {
//            x(i)=xold(i)+alam*p(i);
//        }
//        if(x(0) < 0) x(0) = alpha_start;
//        if(x(1) < 0) x(1) = beta_start;
//         cout<<"alpha, beta 3: "<< x <<endl;
//        idum -= 1;
//        vec Es = MCImportance(idum, x(0),x(1),min_steps, slater, hamiltonian, corr, allEnergies);
//        *f = Es(0)/(min_steps * nParticles);
//        if (alam < alamin) {
//            for (i=0;i<n;i++) x(i)=xold(i);
//            *check=1;
//            return;
//        } else if (*f <= fold+ALF*alam*slope) return;
//        else {
//            if (alam == 1.0)
//                tmplam = -slope/(2.0*(*f-fold-slope));
//            else {
//                rhs1 = *f-fold-alam*slope;
//                rhs2=f2-fold2-alam2*slope;
//                a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
//                b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
//                if (a == 0.0) tmplam = -slope/(2.0*b);
//                else {
//                    disc=b*b-3.0*a*slope;
//                    if (disc<0.0) cout << "Roundoff problem in lnsrch." << endl;
//                    else tmplam=(-b+sqrt(disc))/(3.0*a);
//                }
//                if (tmplam>0.5*alam)
//                    tmplam=0.5*alam;
//            }
//        }
//        alam2=alam;
//        f2 = *f;
//        fold2=fold;
//        alam=FMAX(tmplam,0.1*alam);
//    }
//}
//#undef ALF
//#undef TOLX

