#include "vmcimportance.h"
#include "lib.h"
#include "WaveFunction.h"
#include "hamiltonian.h"


#include <armadillo>
#include <fstream>
#include <iostream>
#include <mpi.h>


using namespace arma;
using namespace std;

VMCImportance::VMCImportance():
    nDimensions(3),
    charge(2),
    nParticles(2),
    h(0.001),
    h2(1000000),
    idum(-1),
    nCycles(1000000),
    alpha_min(1.6),
    alpha_max(1.7),
    alpha_steps(2),
    beta_min(0.3),
    beta_max(0.4),
    beta_steps(2),
    timestep(0.05),
    D(0.5)


{
}

void VMCImportance::runMonteCarloIntegration(int argc, char *argv[])
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

    vec alphas = zeros(alpha_steps);
    vec betas = zeros(beta_steps);

    mat energies = zeros(alpha_steps,beta_steps);
    mat energySquareds = zeros(alpha_steps,beta_steps);

    int id, np;
    //double eTime,sTime;

    //mpd --ncpus=4 &
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
//            sTime = MPI_Wtime();

    int mpi_steps = alpha_steps/np;
    int remainder = alpha_steps%np;

    //cout << "mpi steps " << mpi_steps << " remainder " << remainder << endl;

    mat pEnergies = zeros(alpha_steps,beta_steps);
    mat pEnergySquareds = zeros(alpha_steps,beta_steps);
    mat qForceOld = zeros(alpha_steps,beta_steps);
    mat qForceNew = zeros(alpha_steps,beta_steps);

    int mpi_start = mpi_steps*id;
    int mpi_stop = mpi_start + mpi_steps;
    if (id == np-1) mpi_stop += remainder;

    //cout << "Id " << id << " k start, stop " << mpi_start << " " << mpi_stop << endl;

    for (int k=mpi_start; k<mpi_stop; k++) {
        alpha = alpha_min + k*alpha_step;
        alphas(k) = alpha;
        for (int l=0; l<beta_steps; l++) {
            beta = beta_min + l*beta_step;
            cout << "k,l,alpha,beta: " << k << " " << l <<" "<< alpha << " " << beta <<endl;
            betas(l) = beta;

            // initial positions
            for(int i = 0; i < nParticles; i++) {
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = gaussianDeviate(&idum)*sqrt(timestep);
                }
            }

            waveFunctionOld = function->waveFunction(rOld, alpha, beta);
            qForceOld = quantumForce(rOld, alpha, beta, waveFunctionOld,function);

            // loop over Monte Carlo cycles
            for(int cycle = 0; cycle < nCycles; cycle++) {

                // New position to test
                for(int i = 0; i < nParticles; i++) {
                    for(int j = 0; j < nDimensions; j++) {
                        rNew(i,j) = rOld(i,j) + gaussianDeviate(&idum)*sqrt(timestep) + qForceOld(i,j)*timestep*D;
                    }

                    //Move only one particle.
                    for (int k=0; k<nParticles; k++) {
                        if(k != 0) {
                            for(int j=0; j<nDimensions; j++) {
                                rNew(k,j) = rOld(k,j);
                            }
                        }
                    }

                    // Recalculate the value of the wave function
                    waveFunctionNew = function->waveFunction(rNew, alpha, beta);
                    qForceNew = quantumForce(rNew, alpha, beta, waveFunctionNew,function);

                    //Greens function
                    double greensFunction = 0;
                    for(int j=0; j<nDimensions; j++) {
                        greensFunction += 0.5*(qForceOld(i,j) + qForceNew(i,j)) * (0.5*D*timestep*(qForceOld(i,j) - qForceNew(i,j) - rNew(i,j) + rOld(i,j)));
                    }
                    greensFunction = exp(greensFunction);

                    ++count_total;

                    // Check for step acceptance (if yes, update position, if no, reset position)
                    if(ran2(&idum) <= greensFunction * (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                        ++accepted_steps;
                        for(int j = 0; j < nDimensions; j++) {
                            rOld(i,j) = rNew(i,j);
                            qForceOld(i,j) = qForceNew(i,j);
                            rowvec r12 = rOld.row(1) - rOld.row(0);
                            average_dist += norm(r12, 2);
                        }
                        waveFunctionOld = waveFunctionNew;
                    }
//                    else {
//                        for(int j = 0; j < nDimensions; j++) {
//                            rNew(i,j) = rOld(i,j);
//                            qForceNew(i,j) = qForceOld(i,j);
//                        }
//                    }
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
            cout << "Energy: " << pEnergies(k,l)*2*13.6 << endl;
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

    if(id == 0) {
        cout << energies*2*13.6 << endl;
    printFile(*file_energies, *file_energySquareds, *file_alpha, energies, energySquareds, alphas, betas);
    }


}


mat VMCImportance::quantumForce(const mat &r, double alpha_, double beta_, double wf, WaveFunction *function) {

    mat qforce = zeros(nParticles, nDimensions);
    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    //First derivative

    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = function->waveFunction(rMinus, alpha_, beta_);
            waveFunctionPlus = function->waveFunction(rPlus, alpha_, beta_);
            qforce(i,j) = (waveFunctionPlus - waveFunctionMinus)*2/(wf*2*h);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }

    return qforce;
}

double VMCImportance::gaussianDeviate(long *idum) {

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



void VMCImportance::printFile(const char &file_energies, const char &file_energySquareds, const char &file_alpha, const mat &energies, const mat &energiesSquared, const vec alphas, const vec betas)
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

