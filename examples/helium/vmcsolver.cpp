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
    charge(2),      //Charge of atomic nucleus
    nParticles(2),  //No of electrons in atom
    h(0.001),       //Constants used in numeric derivatives
    h2(1000000),
    idum(-1),       //Random generator seed
    nCycles(1000000),  //No of MC cycles
    alpha_min(1.8), //Alpha, minimum value (loop over alpha values)
    alpha_max(1.8), //Alpha, max value
    alpha_steps(1), //No of steps in alpha loop
    beta_min(0.6),  //Beta, min value
    beta_max(0.6),  //Beta, max value
    beta_steps(1),  //No of steps in beta loop
    timestep(0.01), //Timestep in importance sampling
    D(0.5),         //Constant in importance sampling
    stepLength(1),  //Steplength in brute force Monte Carlo
    minimise_var(false), //Use optimizer to find best values for alpha and beta
    min_steps(50000)//Number of MC cycles for optimizer to run


{
}

//Start the MC method
void VMCSolver::runMonteCarloIntegration(int argc, char *argv[])
{
    bool print_blockdata = false; //Print data to file for blocking (standard deviation)
    char file_energies[] = "../../../output/energy.txt"; //Print data to file (alpha/beta loops)
    char file_energySquareds[] = "../../../output/squareds.txt";
    char file_alpha[] = "../../../output/alpha_beta.txt";

    //Make WaveFunction and Hamiltonian-objects
    WaveFunction *function = new WaveFunction(nParticles, nDimensions);
    Hamiltonian *hamiltonian = new Hamiltonian(nParticles, nDimensions, h, h2, charge);

    double energySum = 0;
    double energySquaredSum = 0;

    double alpha = 0;
    double beta = 0;

    //Compute the steps sizes for the alpha and beta loops
    double alpha_step = (alpha_max - alpha_min)/(alpha_steps-1);
    double beta_step = (beta_max - beta_min)/(beta_steps-1);
    if(alpha_max == alpha_min)  alpha_step = 1;
    if(beta_max == beta_min) beta_step = 1;

    vec alphas = zeros(alpha_steps);
    vec betas = zeros(beta_steps);

    mat energies = zeros(alpha_steps,beta_steps);
    mat energySquareds = zeros(alpha_steps,beta_steps);

    int id, np;

    //Start parallel threads
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    double myTime,mintime, maxtime,avgtime;
    myTime = MPI_Wtime(); //Timing of threads

    int mpi_steps = nCycles/np; //Numver of MC cycles per thread
    idum = idum-id; //Different seed for each thread

    double* allEnergies = new double[mpi_steps+1];

    //Matrices to store results
    mat pEnergies = zeros(alpha_steps,beta_steps);
    mat pEnergySquareds = zeros(alpha_steps,beta_steps);


    for (int k=0; k<alpha_steps; k++) { //Loop over alpha values
        alpha = alpha_min + k*alpha_step;
        alphas(k) = alpha;
        for (int l=0; l<beta_steps; l++) { //Loop over beta values
            beta = beta_min + l*beta_step;
            betas(l) = beta;

            cout << "ID, k,l,alpha,beta: " << id << " "<< k << " " << l <<" "<< alpha << " " << beta <<endl;

            if(minimise_var) { //If optimization of alpha and beta

                double gtol = 1e-4;
                int iter;
                double fret;
                vec p = zeros<vec>(2,1);
                p(0) = alpha;
                p(1) = beta;
                int n = 2;

                vec ans = steepest_descent(idum, p, n, gtol, min_steps, &iter, &fret, hamiltonian, function);
                double alpha_new = ans(0);
                double beta_new = ans(1);

                MPI_Barrier(MPI_COMM_WORLD);

                //Find average values of alpha and beta over threads
                MPI_Allreduce(&alpha_new, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&beta_new, &beta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                alpha = alpha/np;
                beta = beta/np;
                cout << "Final alpha, beta: "<< alpha<<" "<<beta<<endl;

            }

            //Importance sampling or brute force MC:
            vec energyVec = MCImportance(idum, alpha, beta, mpi_steps, function, hamiltonian, allEnergies);
            energySum = energyVec(0);
            energySquaredSum = energyVec(1);
            //MCSampling(alpha, beta, mpi_steps, function, hamiltonian, energySum, energySquaredSum, allEnergies);


            if(print_blockdata) {
                ostringstream ost;
                //ost << "/mn/korona/rp-s1/alborg/4411/helium/examples/vmc-simple/DATA/data" << id << ".mat" ;
                ost << "/home/anette/helium/examples/helium/DATA/data" << id << ".mat" ;

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
            pEnergies(k,l) = energySum/(nCycles * nParticles);
            pEnergySquareds(k,l) = energySquaredSum/(nCycles * nParticles);

            cout << "--------------------------" << endl;

            energySum = 0;
            energySquaredSum = 0;


        } //End beta loop
    } //End alpha loop


    MPI_Barrier(MPI_COMM_WORLD);

    //Gather data from all threads:
    MPI_Allreduce(pEnergies.memptr(), energies.memptr(), alpha_steps*beta_steps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(pEnergySquareds.memptr(), energySquareds.memptr(), alpha_steps*beta_steps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    myTime = MPI_Wtime() - myTime;
    MPI_Reduce(&myTime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&myTime, &mintime, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
    MPI_Reduce(&myTime, &avgtime, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
    MPI_Finalize();
    //End of parallel threads

    if (id == 0) { //Thread no 0:
        cout << "Energies: "<<energies << endl; //*2*13.6
        cout<< "Squareds: "<< energySquareds<<endl;
        cout << "Var: "<<  (energySquareds - square(energies))/nCycles << endl;
        cout << "Std: "<<  sqrt((energySquareds - square(energies))/nCycles) << endl;
        //Print data to file:
       // printFile(*file_energies, *file_energySquareds, *file_alpha, energies, energySquareds, alphas, betas);
        avgtime /= np;
        cout << "Min time: " << mintime << ", max time: " << maxtime << ", avg time: " << avgtime << endl;
    }

    delete[] allEnergies;
}

//Do importance sampling MC:
vec VMCSolver::MCImportance(long idum, double alpha, double beta, int mpi_steps, WaveFunction *function, Hamiltonian *hamiltonian, double *allEnergies) {

    mat qForceOld = zeros(alpha_steps,beta_steps);
    mat qForceNew = zeros(alpha_steps,beta_steps);
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);
    int accepted_steps = 0;
    int count_total = 0;
    double deltaE = 0;
    double waveFunctionOld = 0;
    double waveFunctionNew = 0;
    double cycleE = 0;

    vec deltaPsi = zeros<vec>(2);
    vec deltaPsiE = zeros<vec>(2);
    vec energySums = zeros<vec>(6);

    //Find initial positions for all electrons
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = gaussianDeviate(&idum)*sqrt(timestep);
        }
    }

    //Compute the wavefunction and quantum force
    waveFunctionOld = function->waveFunction(rOld, alpha, beta);
    qForceOld = quantumForce(rOld, alpha, beta, waveFunctionOld,function);

    for(int cycle = 0; cycle < mpi_steps; cycle++) { // loop over Monte Carlo cycles

        for(int i = 0; i < nParticles; i++) { //Loop over particles (electrons)

             // New position for current electron:
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + gaussianDeviate(&idum)*sqrt(timestep) + qForceOld(i,j)*timestep*D;
            }

            //Move only one particle.
            for (int g=0; g<nParticles; g++) {
                if(g != i) {
                    for(int j=0; j<nDimensions; j++) {
                        rNew(g,j) = rOld(g,j);
                    }
                }
            }

            // Recalculate the wave function and quantum force
            waveFunctionNew = function->waveFunction(rNew, alpha, beta);
            qForceNew = quantumForce(rNew, alpha, beta, waveFunctionNew,function);

            //Greens function
            double greensFunction = 0;
            for(int j=0; j<nDimensions; j++) {
                greensFunction += 0.5*(qForceOld(i,j) + qForceNew(i,j)) * (0.5*D*timestep*(qForceOld(i,j) - qForceNew(i,j)) - rNew(i,j) + rOld(i,j));
            }
            greensFunction = exp(greensFunction);

            ++count_total;

            // Check for step acceptance (if yes, update position, if no, reset position)
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

            //Get contribution to energy
            deltaE = hamiltonian->localEnergy(rNew, alpha, beta, function);
            //deltaE = hamiltonian->analyticEnergyHe(rNew, alpha, beta); //Get analytic energy of He
            energySums(0) += deltaE;
            energySums(1) += deltaE*deltaE;
            cycleE +=deltaE;
            if(minimise_var) { //If optimizer is in use, get expectance value of dPsi/dalpha and dPsi/dbeta
                deltaPsi = hamiltonian->dPsi(rNew,alpha,beta,function);
                deltaPsiE(0) = deltaE*deltaPsi(0);
                deltaPsiE(1) = deltaE*deltaPsi(1);
                energySums(2) += deltaPsi(0);
                energySums(3) += deltaPsi(1);
                energySums(4) += deltaPsiE(0);
                energySums(5) += deltaPsiE(1);
            }

        } //End particle loop

        allEnergies[cycle] = cycleE/nParticles; //Store energy for this MC cycle (for blocking method)
        cycleE = 0;

    } //End Monte Carlo loop

    cout << "accepted steps: " << 100*accepted_steps/count_total <<"%"<< endl;

     return energySums;

}

//Do brute force MC
void VMCSolver::MCSampling(double alpha, double beta, int mpi_steps, WaveFunction *function, Hamiltonian *hamiltonian, double &energySum, double &energySquaredSum, double *allEnergies) {

    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);
    int accepted_steps = 0;
    int count_total = 0;
    double deltaE = 0;
    double waveFunctionOld = 0;
    double waveFunctionNew = 0;
    vec average_dist = zeros(3);
    double cycleE = 0;


    //Adjust step length, depending on alpha:
    if(alpha==1) stepLength = 2.5;
    if(alpha>1) stepLength = 2.25;
    if(alpha>=1.2) stepLength = 2;
    if(alpha>=1.4) stepLength = 1.75;
    if(alpha>=1.6) stepLength = 1.5;
    if(alpha>=2) stepLength = 1.25;
    if(alpha>=2.3) stepLength = 1;

    //Get initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
        }
    }
    rNew = rOld;

    for(int cycle = 0; cycle < mpi_steps; cycle++) { // loop over Monte Carlo cycles

        // Store the current value of the wave function
        waveFunctionOld = function->waveFunction(rOld, alpha, beta);


        for(int i = 0; i < nParticles; i++) { //Loop over particles

            // New position for current particle:
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }

            // Recalculate the value of the wave function
            waveFunctionNew = function->waveFunction(rNew, alpha, beta);
            ++count_total;

            // Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                ++accepted_steps;
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                    waveFunctionOld = waveFunctionNew;
                    rowvec r12 = rOld.row(1) - rOld.row(0);
                    double r12norm = 0;
                    double r1norm = 0;
                    double r2norm = 0;
                    for(int d= 0;d<nDimensions;d++) {
                        r12norm += pow(r12(d),2);
                        r1norm += pow(rOld(0,d),2);
                        r2norm += pow(rOld(1,d),2);
                    }
                    average_dist(0) += sqrt(r12norm);
                    average_dist(1) += sqrt(r1norm);
                    average_dist(2) += sqrt(r2norm);
                }
            } else {
                for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j);
                }
            }
            // Get contribution to energy
            deltaE = hamiltonian->localEnergy(rNew, alpha, beta, function);
            cycleE +=deltaE;
            //deltaE = hamiltonian->analyticEnergyHe(rNew, alpha, beta); //Analytic energy He
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;


        } //End particle loop
        allEnergies[cycle] = cycleE/nParticles; //Store energy for this MC cycle (for blocking method)
        cycleE = 0;

    } //End MC loop

    cout << "average distance r1: "<<average_dist(1)/(mpi_steps*nParticles) <<endl;
    cout << "average distance r2: "<<average_dist(2)/(mpi_steps*nParticles) <<endl;
    cout << "average distance r12: "<<average_dist(0)/(mpi_steps*nParticles) <<endl;
    cout << "accepted steps: " << 100*accepted_steps/count_total <<"%"<< endl;
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
    double alpha_step = 1;
    double beta_step = 1;
    int i = 0;
    double test = 0;
    double step_reduce = 2;

    //Get E for current alpha and beta, do MC sample
    vec Es = MCImportance(idum, alpha, beta, min_steps,function,hamiltonian,allEnergies);
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
        cout<<"dE alpha, step: "<<dE(0)<<" "<<alpha_step<<endl;
        dEold = dE;

        //Get E for new alpha and current beta, do MC sample
        Es = MCImportance(idum, alpha_new,beta_new,min_steps, function,hamiltonian, allEnergies);
        Enew = Es(0)/(min_steps * nParticles);
        dPsi(0) = Es(2)/(min_steps * nParticles);
        dPsi(1) = Es(3)/(min_steps * nParticles);
        dPsi_Elocal(0) = Es(4)/(min_steps * nParticles);
        dPsi_Elocal(1) = Es(5)/(min_steps * nParticles);
        dE = gradE(dPsi, E, dPsi_Elocal); //Get derivatives of E wrt alpha and beta
        //If derivatives have changed sign, reduce step length
        if(dE(0)*dEold(0) < 0) alpha_step = alpha_step/step_reduce;
        if(dE(1)*dEold(1) < 0) beta_step = beta_step/step_reduce;

        beta_new = beta - beta_step*dE(1); //Get new value of beta
        if(beta_new < 0) { //If the new beta is negative,
            while(beta_new < 0) { //Reduce step length until new beta is positive
                beta_step = beta_step/step_reduce;
                beta_new = beta - beta_step*dE(1);
            }
        }
        cout<<"dE beta, step: "<<dE(1)<<beta_step<<endl;
        dEold = dE;

        //Get E for new alpha and beta, do MC sample
        Es = MCImportance(idum, alpha_new,beta_new,min_steps, function,hamiltonian, allEnergies);
        Enew = Es(0)/(min_steps * nParticles);
        dPsi(0) = Es(2)/(min_steps * nParticles);
        dPsi(1) = Es(3)/(min_steps * nParticles);
        dPsi_Elocal(0) = Es(4)/(min_steps * nParticles);
        dPsi_Elocal(1) = Es(5)/(min_steps * nParticles);
        dE = gradE(dPsi, E, dPsi_Elocal);//Get derivatives of E wrt alpha and beta
        //If derivatives have changed sign, reduce step length
        if(dE(0)*dEold(0) < 0) alpha_step = alpha_step/step_reduce;
        if(dE(1)*dEold(1) < 0) beta_step = beta_step/step_reduce;

        cout<<"alpha new: "<< alpha_new<<endl;
        cout <<"beta new: "<< beta_new<<" "<<endl;
        cout<<"Enew: "<<Enew<<endl;
        cout <<"----------"<<endl;

        test = abs(Enew-E);
        if(test < gtol) break;  //If change in energy is smaller than tolerance, break out of loop
        E = Enew; //Else: Update E, alpha and beta
        beta = beta_new;
        alpha = alpha_new;

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

//Get quantum force (first derivative of wave function)
mat VMCSolver::quantumForce(const mat &r, double alpha_, double beta_, double wf, WaveFunction *function) {

    mat qforce = zeros(nParticles, nDimensions);
    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    //First derivative, numerical approach

    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) = r(i,j)+h;
            rMinus(i,j) = r(i,j)-h;
            waveFunctionMinus = function->waveFunction(rMinus, alpha_, beta_);
            waveFunctionPlus = function->waveFunction(rPlus, alpha_, beta_);
            qforce(i,j) = (waveFunctionPlus - waveFunctionMinus)/(wf*h);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }

    return qforce;
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


//Print energy data to file
void VMCSolver::printFile(const char &file_energies, const char &file_energySquareds, const char &file_alpha, const mat &energies, const mat &energiesSquared, const vec alphas, const vec betas)
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


}

