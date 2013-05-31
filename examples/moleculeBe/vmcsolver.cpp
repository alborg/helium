#include "vmcsolver.h"
#include "lib.h"
#include "WaveFunction.h"
#include "hamiltonian.h"
#include "slaterdeterminant.h"
#include "correlation.h"


#include <armadillo>
#include <fstream>
#include <iostream>
#include <mpi.h>


using namespace arma;
using namespace std;

VMCSolver::VMCSolver():

  nDimensions(3), //No of dimensions (1D, 2D, 3D, ...)
  h(0.001),       //Constants used in numeric derivatives
  h2(1000000),
  nCycles(1000000),  //No of MC cycles
  timestep(0.01), //Timestep in importance sampling
  D(0.5),         //Constant in importance sampling
  stepLength(1),  //Steplength in brute force Monte Carlo
  minimise_var(false), //Use optimizer to find best values for alpha and beta
  min_steps(50000),//Number of MC cycles for optimizer to run
  alpha(3.7),
  beta(0.23),
  R(1.4),         //Distance between protons (nucleii)
  nProtons(2),    //Total number of nuclei
  nElectrons(1),  //No of electron per atom
  nParticles(nProtons*nElectrons) //Total number of electrons

{
}

void VMCSolver::runMonteCarloIntegration(int argc, char *argv[])
{
    bool printToFile = false; //Print to blocking files

    rProtons = zeros<mat>(nProtons, nDimensions);
    rProtons(0,2) = -R/2;
    rProtons(1,2) = R/2;

    slaterDeterminant *slater = new slaterDeterminant(nDimensions,nProtons,nElectrons);
    Hamiltonian *hamiltonian = new Hamiltonian(nDimensions,h,h2,nProtons,nElectrons,R,rProtons);
    correlation *corr = new correlation(nDimensions,nProtons,nElectrons);

    double energies = 0;
    double energySquareds = 0;
    long idum = -1;
    int id, np;


    //Start parallel threads
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    //Star timing
    double myTime,mintime, maxtime,avgtime;
    myTime = MPI_Wtime();

    int mpi_steps = nCycles/np;  //No of MC cycles per thread
    idum = idum-id;  //Random seed for each thread

    double* allEnergies = new double[mpi_steps+1];
    double pEnergies = 0;
    double pEnergySquareds = 0;
    vec energySums = zeros<vec>(6);

    cout << "ID: " << id << endl;


    if(minimise_var) { //Optimize alpha, beta

        double gtol = 5e-4;
        int iter;
        double fret;
        vec p = zeros<vec>(2,1);
        p(0) = alpha;
        p(1) = beta;
        int n = 2;

        vec ans = steepest_descent(idum, p, n, gtol, min_steps, &iter, &fret, slater,corr, hamiltonian);
        cout <<ans<<endl;
        double alpha_new = ans(0);
        double beta_new = ans(1);

        MPI_Barrier(MPI_COMM_WORLD);
        //Gather new alpha, beta
        MPI_Allreduce(&alpha_new, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&beta_new, &beta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        alpha = alpha/np;
        beta = beta/np;
        cout << "Final alpha, beta: "<< alpha<<" "<<beta<<endl;


    }

    //Call importance sampling MC
    energySums = MCImportance(idum, alpha, beta, mpi_steps, slater, hamiltonian, corr, allEnergies);

    if(printToFile) { //Blocking files
        ostringstream ost;
        ost << "/DATA/data" << id << ".mat" ;

        ofstream blockofile;
        blockofile.open( ost.str( ).c_str( ),ios::out | ios::binary );
        if (blockofile.is_open())
        {
            blockofile.write((char*)(allEnergies+1) , mpi_steps*sizeof(double)) ;
            blockofile.close();
        }
        else cout << "Unable to open data file for process " << id << endl;
    }

    //Find average values of energies:
    pEnergies = energySums(0)/(nCycles * nParticles);
    pEnergySquareds = energySums(1)/(nCycles * nParticles);

    cout << "--------------------------" << endl;



    MPI_Barrier(MPI_COMM_WORLD);

    //Gather energy data from threads
    MPI_Allreduce(&pEnergies, &energies, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&pEnergySquareds, &energySquareds, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    myTime = MPI_Wtime() - myTime;
    MPI_Reduce(&myTime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&myTime, &mintime, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
    MPI_Reduce(&myTime, &avgtime, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
    MPI_Finalize();
    //Parallel threads ended

    if (id == 0) {
        cout << "Energy: " << energies << endl; //*2*13.6
        cout << "Energy squareds: " << energySquareds << endl; //*2*13.6*2*13.6
        avgtime /= np;
        cout << "Min time: " << mintime << ", max time: " << maxtime << ", avg time: " << avgtime << endl;
    }

    delete[] allEnergies;
}

//Importance sampling MC
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


    //Get initial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = gaussianDeviate(&idum)*sqrt(timestep);
        }
    }


    //Build the full Slater matrix and invert (done only once):
    slater->buildDeterminant(rOld, alpha);

    for(int cycle = 0; cycle < mpi_steps; cycle++) { // loop over Monte Carlo cycles

        for(int i = 0; i < nParticles; i++) { //Loop over particles

            ratio = 1.0;
            qForceOld = slater->gradientWaveFunction(rOld, i, ratio, alpha); //Get quantum force

            // New position current particle
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

            //Get the ratio of the new to the old determinant (wavefunction), Slater and correlation factor.
            ratio = slater->getRatioDeterminant(i, rNew, alpha);
            ratioCorr = corr->getRatioJastrow(i, rOld, rNew, beta);

            //Get new quantum force
            qForceNew = slater->gradientWaveFunction(rNew, i, ratio, alpha);

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
               slater->updateDeterminant(rNew, rOld, i, alpha, ratio); //Update the Slater determinant
               for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);   
                }
            }
            else {
               for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j);
                }
            }


           //Get contribution to energy
            deltaE = hamiltonian->localEnergy(rNew, alpha, beta, slater,corr);
            energySums(0) += deltaE;
            energySums(1) += deltaE*deltaE;
            allEnergies[cycle] += deltaE;
            cycleE += deltaE;
            if(minimise_var) { //If optimizer is in use, get expectance value of dPsi/dalpha and dPsi/dbeta
                deltaPsi = hamiltonian->dPsi(rNew,alpha,beta,slater,corr);
                deltaPsiE(0) = deltaE*deltaPsi(0);
                deltaPsiE(1) = deltaE*deltaPsi(1);
                energySums(2) += deltaPsi(0);
                energySums(3) += deltaPsi(1);
                energySums(4) += deltaPsiE(0);
                energySums(5) += deltaPsiE(1);
            }


        } //End particle loop

        allEnergies[cycle] += cycleE; //Store energy for this MC cycle (for blocking method)
        cycleE = 0;

    } //End Monte Carlo loop

    cout << "accepted steps: " << 100*accepted_steps/count_total << "%" << endl;

    return energySums;

}



//Get random numbers with a Gaussian pdf
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


//Optimization of alpha and beta:
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
    double alpha_step = 1;
    double beta_step = 1;
    int i = 0;
    double test;
    double step_reduce = 2;

    //Get E for current alpha and beta, do MC sample
    vec Es = MCImportance(idum, alpha,beta,min_steps, slater, hamiltonian, corr, allEnergies);
    E = Es(0)/(min_steps * nParticles);
    dPsi(0) = Es(2)/(min_steps * nParticles);
    dPsi(1) = Es(3)/(min_steps * nParticles);
    dPsi_Elocal(0) = Es(4)/(min_steps * nParticles);
    dPsi_Elocal(1) = Es(5)/(min_steps * nParticles);
    dE = gradE(dPsi, E, dPsi_Elocal);  //Get derivatives of E wrt alpha and beta

    cout <<"E: "<<E<<endl;

    while(i<maxIter) { //Loop until enough iterations

        alpha_new = alpha - alpha_step*dE(0); //Get new value of alpha
        if(alpha_new < 0) { //If the new alpha is negative,
            while(alpha_new < 0) { //Reduce step length until new alpha is positive
                alpha_step = alpha_step/step_reduce;
                alpha_new = alpha - alpha_step*dE(0);
            }
        }
        cout<<"dE alpha: "<<dE(0)<<endl;
        dEold = dE;

        //Get E for current alpha and beta, do MC sample
        Es = MCImportance(idum, alpha_new,beta,min_steps, slater, hamiltonian, corr, allEnergies);
        Enew = Es(0)/(min_steps * nParticles);
        dPsi(0) = Es(2)/(min_steps * nParticles);
        dPsi(1) = Es(3)/(min_steps * nParticles);
        dPsi_Elocal(0) = Es(4)/(min_steps * nParticles);
        dPsi_Elocal(1) = Es(5)/(min_steps * nParticles);
        dE = gradE(dPsi, E, dPsi_Elocal);  //Get derivatives of E wrt alpha and beta
        //If derivatives have changed sign, reduce step length
        if(dE(0)*dEold(0) < 0) alpha_step = alpha_step/step_reduce;
        if(dE(1)*dEold(1) < 0) beta_step = beta_step/step_reduce;
        cout <<"Enew: "<<Enew<<endl;

        test = abs(Enew-E);
        if(test < gtol) break;
        E = Enew;

        cout <<"Alpha new: "<< alpha_new <<endl;
        cout <<"dE, Step: "<< dEold(0)<<" "<<alpha_step << endl;
        cout<<"Enew: "<<Enew<<endl;

        beta_new = beta - beta_step*dE(1); //Get new value of beta
        if(beta_new < 0) { //If the new beta is negative,
            while(beta_new < 0) { //Reduce step length until new beta is positive
                beta_step = beta_step/step_reduce;
                beta_new = beta - beta_step*dE(1);
            }
        }
        dEold = dE;

        //Get E for current alpha and beta, do MC sample
        Es = MCImportance(idum, alpha_new,beta_new,min_steps, slater, hamiltonian, corr, allEnergies);
        Enew = Es(0)/(min_steps * nParticles);
        dPsi(0) = Es(2)/(min_steps * nParticles);
        dPsi(1) = Es(3)/(min_steps * nParticles);
        dPsi_Elocal(0) = Es(4)/(min_steps * nParticles);
        dPsi_Elocal(1) = Es(5)/(min_steps * nParticles);
        dE = gradE(dPsi, E, dPsi_Elocal);  //Get derivatives of E wrt alpha and beta
        //If derivatives have changed sign, reduce step length
        if(dE(0)*dEold(0) < 0) alpha_step = alpha_step/step_reduce;
        if(dE(1)*dEold(1) < 0) beta_step = beta_step/step_reduce;

        cout <<"beta new: "<< beta_new<<" "<<endl;
        cout <<"dE, Step: "<< dEold(1)<<" "<<beta_step<<" "<<endl;
        cout<<"Enew: "<<Enew<<endl;
        cout <<"----------"<<endl;

        test = abs(Enew-E);
        if(test < gtol) break;  //If change in energy is smaller than tolerance, break out of loop
        E = Enew; //Else: Update E, alpha and beta
        alpha = alpha_new;
        beta = beta_new;

        i++;
    }

    answers(0) = alpha_new;
    answers(1) = beta_new;
    answers(2) = Enew;
    answers(3) = i;

    return answers;
}

//Get derivatives of energy E wrt alpha, beta
vec VMCSolver::gradE(vec dPsi, double Elocal, vec dPsi_Elocal) {

    vec dE = zeros<vec>(2);
    dE(0) = 2*(dPsi_Elocal(0) - dPsi(0)*Elocal);
    dE(1) = 2*(dPsi_Elocal(1) - dPsi(1)*Elocal);

    return dE;
}

