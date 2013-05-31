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

    nDimensions(3), //No of dimensions (1D, 2D, 3D, ...)
    charge(1),      //Charge of atomic nucleus
    h(0.001),       //Constants used in numeric derivatives
    h2(1000000),
    nCycles(1000000),  //No of MC cycles
    timestep(0.01), //Timestep in importance sampling
    D(0.5),         //Constant in importance sampling
    stepLength(1),  //Steplength in brute force Monte Carlo
    minimise_var(false), //Use optimizer to find best values for alpha and beta
    min_steps(50000),//Number of MC cycles for optimizer to run
    alpha(1.3),
    beta(0.3),
    R(1.4),         //Distance between protons (nucleii)
    nProtons(2),    //Total number of nuclei
    nElectrons(1),  //No of electron per atom
    nParticles(nProtons*nElectrons), //Total number of electrons
    printToFile(false) //Blocking

{
}

void VMCSolver::runMonteCarloIntegration(int argc, char *argv[])
{

    WaveFunction *function = new WaveFunction(nDimensions,nProtons,nElectrons);
    Hamiltonian *hamiltonian = new Hamiltonian(nProtons, nElectrons,nDimensions,h,h2,charge);

    double energies = 0;
    double energySquareds = 0;
    long idum = -1;
    mat Rs = zeros<vec>(20,4);

    int id, np;

    rProtons = zeros<mat>(nProtons, nDimensions);

    //Start parallel threads
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    //Start timing
    double myTime,mintime, maxtime,avgtime;
    myTime = MPI_Wtime();

    //No of MC cycles per thread
    int mpi_steps = nCycles/np;
    idum = idum-id; //Random seed for each thread

    double* allEnergies = new double[mpi_steps+1];
    double pEnergies = 0;
    double pEnergySquareds = 0;
    vec energySums = zeros<vec>(2);

    cout << "ID: " << id << endl;

    for(int b=0;b<1;b++) { //Loop over R values, invalidated

        R += 0.2;
        cout<<"R: "<<R<<endl;
        //Set position of protons
        rProtons(0,2) = -R/2;
        rProtons(1,2) = R/2;

        double alpha_new = alpha;
        double beta_new = beta;

        if(minimise_var) { //Optimize alpha and beta

            double gtol = 1e-4;
            int iter;
            double fret;
            vec p = zeros<vec>(2,1);
            p(0) = alpha;
            p(1) = beta;
            int n = 2;

            vec ans = steepest_descent(idum, p, n, gtol, min_steps, &iter, &fret, hamiltonian,function);
            cout <<ans<<endl;
            alpha_new = ans(0);
            beta_new = ans(1);

            MPI_Barrier(MPI_COMM_WORLD);
            //Gather new alpha, beta
            MPI_Allreduce(&alpha_new, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&beta_new, &beta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            alpha = alpha/np;
            beta = beta/np;
            cout << "Final alpha, beta: "<< alpha<<" "<<beta<<endl;

        }

        //Run importance sampling MC
        energySums = MCImportance(idum, alpha_new, beta_new, mpi_steps, hamiltonian, function, allEnergies);

        if(printToFile) {  //Blocking
            cout<<"Printing to blocking file"<<endl;
            ostringstream ost;
            ost << "/home/anette/helium/examples/molecules/DATA/data" << id << ".mat" ;

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

        Rs(b,0) = R;
        Rs(b,1) = alpha_new;
        Rs(b,2) = beta_new;
        Rs(b,3) = pEnergies;
        cout << "--------------------------" << endl;
    }

    cout<<Rs<<endl;

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
        cout << "Energies: " << energies << endl; //*2*13.6
        cout << "Energy squareds: " << energySquareds << endl; //*2*13.6*2*13.6
        avgtime /= np;
        cout << "Min time: " << mintime << ", max time: " << maxtime << ", avg time: " << avgtime << endl;
    }

    delete[] allEnergies;
}

//Importance sampling MC
vec VMCSolver::MCImportance(long idum, double alpha, double beta, int mpi_steps, Hamiltonian *hamiltonian,
                            WaveFunction *function, double *allEnergies) {


    mat rOld = zeros<mat>(nParticles, nDimensions);
    mat rNew = zeros<mat>(nParticles, nDimensions);
    double accepted_steps = 0;
    double count_total = 0;
    double deltaE = 0;
    vec deltaPsi = zeros<vec>(2);
    vec deltaPsiE = zeros<vec>(2);
    double cycleE = 0;
    vec energySums = zeros<vec>(6);
    double r12 = 0;
    double r1tot = 0;
    double r2tot = 0;

    mat qForceOld = zeros(nParticles,nDimensions);
    mat qForceNew = zeros(nParticles,nDimensions);
    double waveFunctionOld = 0;
    double waveFunctionNew = 0;


    //Get initial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = gaussianDeviate(&idum)*sqrt(timestep);
        }
    }

    //Compute the wavefunction and quantum force
    waveFunctionOld = function->waveFunction(rOld, rProtons, alpha, beta);
    qForceOld = quantumForce(rOld, rProtons, alpha, beta, waveFunctionOld,function);

    for(int cycle = 0; cycle < mpi_steps; cycle++) { // loop over Monte Carlo cycles

        for(int i = 0; i < nParticles; i++) { //Loop over particles

            // New position current particle
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

            // Recalculate the wave function and quantum force
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
               r12 += sqrt(pow(rOld(0,0)-rOld(1,0),2)+pow(rOld(0,1)-rOld(1,1),2)+pow(rOld(0,2)-rOld(1,2),2));
               r1tot += sqrt(pow(rOld(0,0),2) + pow(rOld(0,1),2) + pow(rOld(0,2),2));
               r2tot += sqrt(pow(rOld(1,0),2) + pow(rOld(1,1),2) + pow(rOld(1,2),2));
            }
            else {
               for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j);
                    qForceNew(i,j) = qForceOld(i,j);
                }
            }


            //Get contribution to energy
            deltaE = hamiltonian->localEnergy(R, rNew, rProtons, alpha, beta, function);
            energySums(0) += deltaE;
            energySums(1) += deltaE*deltaE;
            allEnergies[cycle] += deltaE;
            cycleE += deltaE;
            if(minimise_var) { //If optimizer is in use, get expectance value of dPsi/dalpha and dPsi/dbeta
                deltaPsi = hamiltonian->dPsi(rNew,rProtons,alpha,beta,function);
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

    cout << "r1: "<<r1tot/(mpi_steps*nParticles)<<endl;
    cout << "r2: "<<r2tot/(mpi_steps*nParticles)<<endl;
    cout << "r12: "<<r12/(mpi_steps*nParticles)<<endl;
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

//Get quantum force (first derivative of wave function)
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

//Optimization of alpha and beta:
vec VMCSolver::steepest_descent(long idum, vec &p, int n, double gtol, int min_steps, int *iter, double *fret,
                       Hamiltonian *hamiltonian, WaveFunction *function)
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
    int maxIter = 50;
    vec answers = zeros<vec>(n+2);
    double E = 0;
    double Enew = 0;
    double alpha_step = 0.1;
    double beta_step = 1;
    int i = 0;
    int j = 0;
    double test;
    double step_reduce = 2;

    double a = 1 + exp(-R/alpha);
    double delta = a - alpha;

    //Optimize alpha using relation alpha = 1 + exp(-R/alpha)
    while(abs(delta) > gtol && j<maxIter) {
        if(delta > 0) alpha_new = alpha + alpha_step;
        else alpha_new = alpha - alpha_step;
        alpha_step = alpha_step/1.2;
        a = 1 + exp(-R/alpha_new);
        cout<<"Alpha, da: "<<alpha_new<<" "<<a<<endl;
        delta = a - alpha_new;
        alpha = alpha_new;
        j++;
    }

    cout <<"Alpha: "<<alpha_new<<endl;

    //Optimize beta:

    //Get E for current alpha and beta, do MC sample
    vec Es = MCImportance(idum, alpha, beta, min_steps,hamiltonian,function,allEnergies);
    E = Es(0)/(min_steps * nParticles);
    dPsi(0) = Es(2)/(min_steps * nParticles);
    dPsi(1) = Es(3)/(min_steps * nParticles);
    dPsi_Elocal(0) = Es(4)/(min_steps * nParticles);
    dPsi_Elocal(1) = Es(5)/(min_steps * nParticles);
    dE = gradE(dPsi, E, dPsi_Elocal); //Get derivatives of E wrt alpha and beta

    cout <<"E: "<<E<<endl;

    while(i<maxIter) { //Loop until enough iterations

        beta_new = beta - beta_step*dE(1); //Get new value of beta
        if(beta_new < 0) { //If the new beta is negative,
            while(beta_new < 0) { //Reduce step length until new beta is positive
                beta_step = beta_step/step_reduce;
                beta_new = beta - beta_step*dE(1);
            }
        }
        cout<<"dE beta: "<<dE(1)<<endl;
        dEold = dE;

        //Get E for current alpha and beta, do MC sample
        Es = MCImportance(idum, alpha_new,beta_new,min_steps,hamiltonian, function, allEnergies);
        Enew = Es(0)/(min_steps * nParticles);
        dPsi(0) = Es(2)/(min_steps * nParticles);
        dPsi(1) = Es(3)/(min_steps * nParticles);
        dPsi_Elocal(0) = Es(4)/(min_steps * nParticles);
        dPsi_Elocal(1) = Es(5)/(min_steps * nParticles);
        dE = gradE(dPsi, E, dPsi_Elocal); //Get derivatives of E wrt alpha and beta
        //If derivatives have changed sign, reduce step length
        if(dE(0)*dEold(0) < 0) alpha_step = alpha_step/step_reduce;
        if(dE(1)*dEold(1) < 0) beta_step = beta_step/step_reduce;

        cout <<"beta new: "<< beta_new<<" "<<endl;
        cout <<"dE, Step: "<< dEold(1)<<" "<<beta_step<<" "<<endl;
        cout<<"Enew: "<<Enew<<endl;
        cout <<"----------"<<endl;

        test = abs(Enew-E);
        if(test < gtol) break; //If change in energy is smaller than tolerance, break out of loop
        E = Enew; //Else: Update E, alpha and beta
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


