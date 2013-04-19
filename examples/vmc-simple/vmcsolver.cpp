#include "vmcsolver.h"
#include "lib.h"
#include "WaveFunction.h"
#include "hamiltonian.h"
#include "slaterdeterminant.h"


#include <armadillo>
#include <fstream>
#include <iostream>
#include <mpi.h>


using namespace arma;
using namespace std;

VMCSolver::VMCSolver():
    nDimensions(3),
    charge(4),
    nParticles(4),
    h(0.001),
    h2(1000000),
    idum(-1),
    nCycles(1000000),
    alpha_min(4),
    alpha_max(4),
    alpha_steps(1),
    beta_min(0.2),
    beta_max(0.2),
    beta_steps(1),
    timestep(0.05),
    D(0.5),
    stepLength(1.0)


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
    WaveFunction *function = new WaveFunction(nParticles, nDimensions);
    Hamiltonian *hamiltonian = new Hamiltonian(nParticles, nDimensions, h, h2, charge);

    double energySum = 0;
    double energySquaredSum = 0;

    double alpha = 0;
    double beta = 0;

    double alpha_step = (alpha_max - alpha_min)/(alpha_steps-1);
    double beta_step = (beta_max - beta_min)/(beta_steps-1);
    if(alpha_max == alpha_min)  alpha_step = 1;
    if(beta_max == beta_min) beta_step = 1;

    vec alphas = zeros(alpha_steps);
    vec betas = zeros(beta_steps);

    mat energies = zeros(alpha_steps,beta_steps);
    mat energySquareds = zeros(alpha_steps,beta_steps);

    int id, np;

    //mpd --ncpus=4 &
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    double myTime,mintime, maxtime,avgtime;
    myTime = MPI_Wtime();


    int mpi_steps = nCycles/np;
    idum = idum-id*0.1;

    double* allEnergies = new double[mpi_steps+1];

    mat pEnergies = zeros(alpha_steps,beta_steps);
    mat pEnergySquareds = zeros(alpha_steps,beta_steps);


    for (int k=0; k<alpha_steps; k++) {
        alpha = alpha_min + k*alpha_step;
        alphas(k) = alpha;
        for (int l=0; l<beta_steps; l++) {
            beta = beta_min + l*beta_step;
            betas(l) = beta;

            cout << "ID, k,l,alpha,beta: " << id << " "<< k << " " << l <<" "<< alpha << " " << beta <<endl;

            //MCImportance(alpha, beta, mpi_steps, function, slater, hamiltonian, energySum, energySquaredSum, allEnergies);
            MCSampling(alpha, beta, mpi_steps, function, slater, hamiltonian, energySum, energySquaredSum, allEnergies);

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

            pEnergies(k,l) = energySum/(nCycles * nParticles);
            pEnergySquareds(k,l) = energySquaredSum/(nCycles * nParticles);

            cout << "--------------------------" << endl;

            energySum = 0;
            energySquaredSum = 0;


        } //End beta loop
    } //End alpha loop


    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allreduce(pEnergies.memptr(), energies.memptr(), alpha_steps*beta_steps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(pEnergySquareds.memptr(), energySquareds.memptr(), alpha_steps*beta_steps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    myTime = MPI_Wtime() - myTime;
    MPI_Reduce(&myTime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&myTime, &mintime, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
    MPI_Reduce(&myTime, &avgtime, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
    MPI_Finalize();

    if (id == 0) {
        cout << "Energies: " << energies << endl; //*2*13.6
        cout << "Energy squareds: " << energySquareds << endl; //*2*13.6*2*13.6
        printFile(*file_energies, *file_energySquareds, *file_alpha, *file_sigma, energies, energySquareds, alphas, betas);
        avgtime /= np;
        cout << "Min time: " << mintime << ", max time: " << maxtime << ", avg time: " << avgtime << endl;
    }

    delete[] allEnergies;
}

void VMCSolver::MCImportance(double alpha, double beta, int mpi_steps, WaveFunction *function, slaterDeterminant *slater, Hamiltonian *hamiltonian, double &energySum, double &energySquaredSum, double *allEnergies) {

    mat qForceOld = zeros(alpha_steps,beta_steps);
    mat qForceNew = zeros(alpha_steps,beta_steps);
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);
    double accepted_steps = 0;
    double count_total = 0;
    double deltaE = 0;
    double waveFunctionOld = 0;
    double waveFunctionNew = 0;

    // initial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = gaussianDeviate(&idum)*sqrt(timestep);
        }
    }

    slater->buildDeterminant(rOld, alpha, beta);
    //waveFunctionOld = function->waveFunction(rOld, alpha, beta);
    waveFunctionOld = slater->getDeterminant();
    qForceOld = quantumForce(rOld, alpha, beta, waveFunctionOld,function);

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < mpi_steps; cycle++) {


        for(int i = 0; i < nParticles; i++) { //Particle

            for(int j = 0; j < nDimensions; j++) {
                // New position to test
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

            // Recalculate the wave function

            //waveFunctionNew = function->waveFunction(rNew, alpha, beta);


            //qForceNew = quantumForce(rNew, alpha, beta, waveFunctionNew,function);

            //Greens function
            //double greensFunction = 0;
            //for(int j=0; j<nDimensions; j++) {
            //    greensFunction += 0.5*(qForceOld(i,j) + qForceNew(i,j)) * (0.5*D*timestep*(qForceOld(i,j) - qForceNew(i,j)) - rNew(i,j) + rOld(i,j));
            //}
            //greensFunction = exp(greensFunction);

            ++count_total;

            // Check for step acceptance (if yes, update position, if no, reset position)
            //if(ran2(&idum) <= greensFunction * (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {

            double ratio = slater->getRatioDeterminant(i, rNew, alpha, beta);

            if(ran2(&idum) <= ratio) {
                ++accepted_steps;
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                    //qForceOld(i,j) = qForceNew(i,j);
                    slater->updateDeterminant(rNew, rOld, i, alpha, beta, ratio);
                    waveFunctionOld = slater->getDeterminant();
                    qForceOld = quantumForce(rOld, alpha, beta, waveFunctionOld,function);

                }
                waveFunctionOld = waveFunctionNew;
            }
            else {
                for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j);
                    //qForceNew(i,j) = qForceOld(i,j);
                }
            }


            // update energies
            deltaE = hamiltonian->localEnergy(rNew, alpha, beta, slater);
            //deltaE = hamiltonian->analyticEnergyH(rNew, alpha, beta);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
            allEnergies[cycle] = deltaE;

        } //End particle loop

    } //End Monte Carlo loop

    cout << "accepted steps: " << 100*accepted_steps/count_total << "%" << endl;

}


void VMCSolver::MCSampling(double alpha, double beta, int mpi_steps, WaveFunction *function, slaterDeterminant *slater, Hamiltonian *hamiltonian, double &energySum, double &energySquaredSum, double *allEnergies) {

    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);
    int accepted_steps = 0;
    int count_total = 0;
    double deltaE = 0;

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

            double ratio = slater->getRatioDeterminant(i, rNew, alpha, beta);


            // Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= ratio) {
                ++accepted_steps;
                for(int d = 0; d < nDimensions; d++) {
                    rOld(i,d) = rNew(i,d);
                }
                slater->updateDeterminant(rNew, rOld, i, alpha, beta, ratio);
            }
            else {
                for(int d = 0; d < nDimensions; d++) {
                    rNew(i,d) = rOld(i,d);
                }
            }

            // update energies
            deltaE = hamiltonian->localEnergy(rNew, alpha, beta, slater);
            //deltaE = hamiltonian->analyticEnergyH(rNew, alpha, beta);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
            allEnergies[cycle] = deltaE;

        } //Particles

    } //Monte Carlo cycles

    cout << "accepted steps: " << 100*accepted_steps/count_total << "%" << endl;
}


mat VMCSolver::quantumForce(const mat &r, double alpha_, double beta_, double wf, WaveFunction *function) {

    mat qforce = zeros(nParticles, nDimensions);
    vec del = zeros<vec>(3,1);

    for(int i = 0; i < nParticles; i++) {

        del = function->gradientWaveFunction(r, i, alpha_, beta_);

        for(int j = 0; j < nDimensions; j++) {

            qforce(i,j) = del(j)/wf;

        }
    }

    return qforce;
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

