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

VMCSolver::VMCSolver() :
    nDimensions(3),
    charge(2),
    stepLength(1.0),
    nParticles(2),
    h(0.001),
    h2(1000000),
    idum(-1),
    nCycles(100000000),
    alpha_min(1.8),
    alpha_max(1.8),
    alpha_steps(1),
    beta_min(0.7),
    beta_max(0.7),
    beta_steps(1)


{
}

void VMCSolver::runMonteCarloIntegration(int argc, char *argv[])
{

    char file_energies[] = "../../../output/energy.txt";
    char file_energySquareds[] = "../../../output/squareds.txt";
    char file_alpha[] = "../../../output/alpha_beta.txt";

    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    WaveFunction *function = new WaveFunction(nParticles, nDimensions);
    Hamiltonian *hamiltonian = new Hamiltonian(nParticles, nDimensions, h, h2, charge);

    double waveFunctionOld = 0;
    double waveFunctionNew = 0;

    double energySum = 0;
    double energySquaredSum = 0;

    int accepted_steps = 0;
    int count_total = 0;

    double deltaE;
    double average_dist = 0;

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
    //double eTime,sTime;

    //mpd --ncpus=4 &
    //mpirun -np 2 exec
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
//            sTime = MPI_Wtime();

    int mpi_steps = nCycles/np;
    idum = idum-id*0.1;

//    int mpi_steps = alpha_steps/np;
//    int remainder = alpha_steps%np;

    //cout << "mpi steps " << mpi_steps << " remainder " << remainder << endl;

    mat pEnergies = zeros(alpha_steps,beta_steps);
    mat pEnergySquareds = zeros(alpha_steps,beta_steps);

//    int mpi_start = mpi_steps*id;
//    int mpi_stop = mpi_start + mpi_steps;
//    if (id == np-1) mpi_stop += remainder;

    //cout << "Id " << id << " k start, stop " << mpi_start << " " << mpi_stop << endl;

     for (int k=0; k<alpha_steps; k++) {
        alpha = alpha_min + k*alpha_step;
        alphas(k) = alpha;
        for (int l=0; l<beta_steps; l++) {
            beta = beta_min + l*beta_step;
            cout << "k,l,alpha,beta: " << k << " " << l <<" "<< alpha << " " << beta <<endl;
            betas(l) = beta;

            // initial trial positions
            for(int i = 0; i < nParticles; i++) {
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
                }
            }
            rNew = rOld;


            // loop over Monte Carlo cycles
            for(int cycle = 0; cycle < mpi_steps; cycle++) {

                // Store the current value of the wave function
                waveFunctionOld = function->waveFunction(rOld, alpha, beta);

                // New position to test
                for(int i = 0; i < nParticles; i++) {
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
                            average_dist += norm(r12, 2);
                        }
                    } else {
                        for(int j = 0; j < nDimensions; j++) {
                            rNew(i,j) = rOld(i,j);
                        }
                    }
                    // update energies
                    deltaE = hamiltonian->localEnergy(rNew, alpha, beta, function);
                    //deltaE = hamiltonian->analyticLocalEnergy(rNew, alpha, beta);
                    energySum += deltaE;
                    energySquaredSum += deltaE*deltaE;

                }

            }


            cout << "Process id: " << id << " " << k + l <<endl;
            cout << "accepted steps, total steps: " << accepted_steps << " " << count_total << endl;

            pEnergies(k,l) = energySum/(nCycles * nParticles);
            pEnergySquareds(k,l) = energySquaredSum/(nCycles * nParticles);
            average_dist = average_dist/accepted_steps;

            cout << "Average r12: " << average_dist << endl;
            //cout << "Energy: " << pEnergies(k,l)*2*13.6 << endl;
            cout << "--------------------------" << endl;

            energySum = 0;
            energySquaredSum = 0;
            accepted_steps = 0;
            count_total = 0;
            average_dist = 0;

        }
    }



//            eTime = MPI_Wtime();
//            pTime = fabs(eTime - sTime);


     MPI_Barrier(MPI_COMM_WORLD);

     MPI_Allreduce(pEnergies.memptr(), energies.memptr(), alpha_steps*beta_steps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
     MPI_Allreduce(pEnergySquareds.memptr(), energySquareds.memptr(), alpha_steps*beta_steps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

     MPI_Finalize();

     if (id == 0) {
         cout << energies*2*13.6 << endl;
         cout << energySquareds << endl;
         cout << energies*energies << endl;
         mat variance = (1.0/nCycles)*(energySquareds - energies*energies);
         mat sigma = sqrt(variance);
         cout <<  sigma << endl;
         printFile(*file_energies, *file_energySquareds, *file_alpha, energies, energySquareds, alphas, betas);
     }

//     MPI_Barrier(MPI_COMM_WORLD);

//    MPI_Allreduce(pEnergies.memptr(), energies.memptr(), alpha_steps*beta_steps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

//    MPI_Allreduce(pEnergySquareds.memptr(), energySquareds.memptr(), alpha_steps*beta_steps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


//    MPI_Finalize();

//    if(id == 0) {
//        cout << energies*2*13.6 << endl;
//    printFile(*file_energies, *file_energySquareds, *file_alpha, energies, energySquareds, alphas, betas);
//    }
}




void VMCSolver::printFile(const char &file_energies, const char &file_energySquareds, const char &file_alpha, const mat &energies, const mat &energiesSquared, const vec alphas, const vec betas)
{

    ofstream myfile(&file_energies);
    if (myfile.is_open())
    {
        for (unsigned int f=0; f<energies.n_rows; f++)
        {
            for (unsigned int l=0; l<energies.n_cols; l++) {
                myfile << energies(f,l)*2*13.6 << " ";
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
                myfile3 << energiesSquared(f,l)*2*13.6 << " ";
            }
            myfile3 << endl;
        }

        myfile3.close();
    }
    else cout << "Unable to open file" << endl;


}
