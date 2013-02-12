#include "vmcsolver.h"
#include "lib.h"

#include <armadillo>
#include <fstream>
#include <iostream>


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
    nCycles(1000000),
    alpha_min(1.4),
    alpha_max(1.9),
    alpha_steps(10),
    beta_min(0.1),
    beta_max(1.1),
    beta_steps(10)


{
}

void VMCSolver::runMonteCarloIntegration()
{

    char file_energies[] = "..\\..\\..\\output\\energy.txt";
    char file_energySquareds[] = "..\\..\\..\\output\\squareds.txt";
    char file_alpha[] = "..\\..\\..\\output\\alpha_beta.txt";

    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    double waveFunctionOld = 0;
    double waveFunctionNew = 0;

    double energySum = 0;
    double energySquaredSum = 0;

    int accepted_steps = 0;
    int count_total = 0;

    double deltaE;

    double alpha = 0;
    double beta = 0;

    double alpha_step = (alpha_max - alpha_min)/(alpha_steps-1);
    double beta_step = (beta_max - beta_min)/(beta_steps-1);

    vec alphas = zeros(alpha_steps);
    vec betas = zeros(beta_steps);

    mat energies = zeros(alpha_steps,beta_steps);
    mat energySquareds = zeros(alpha_steps,beta_steps);


    for (double k=0; k<alpha_steps; k++) {
        alpha = alpha_min + k*alpha_step;
        alphas(k) = alpha;
        for (double l=0; l<beta_steps; l++) {
            beta = beta_min + l*beta_step;
            cout << "k,l,alpha,beta" << k << " " << l <<" "<< alpha << " " << beta <<endl;
            betas(l) = beta;

            // initial trial positions
            for(int i = 0; i < nParticles; i++) {
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
                }
            }
            rNew = rOld;

            // loop over Monte Carlo cycles
            for(int cycle = 0; cycle < nCycles; cycle++) {

                // Store the current value of the wave function
                waveFunctionOld = waveFunction(rOld, alpha, beta);

                // New position to test
                for(int i = 0; i < nParticles; i++) {
                    for(int j = 0; j < nDimensions; j++) {
                        rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
                    }

                    // Recalculate the value of the wave function
                    waveFunctionNew = waveFunction(rNew, alpha, beta);
                    ++count_total;

                    // Check for step acceptance (if yes, update position, if no, reset position)
                    if(ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                        ++accepted_steps;
                        for(int j = 0; j < nDimensions; j++) {
                            rOld(i,j) = rNew(i,j);
                            waveFunctionOld = waveFunctionNew;
                        }
                    } else {
                        for(int j = 0; j < nDimensions; j++) {
                            rNew(i,j) = rOld(i,j);
                        }
                    }
                    // update energies
                    deltaE = localEnergy(rNew, alpha, beta);
                    energySum += deltaE;
                    energySquaredSum += deltaE*deltaE;

                }

            }


            cout << "accepted steps, total steps: " << accepted_steps << " " << count_total << endl;

            energies(k,l) = energySum/(nCycles * nParticles);
            energySquareds(k,l) = energySquaredSum/(nCycles * nParticles);

            cout << "Energy" << energies(k,l) << endl;

            energySum = 0;
            energySquaredSum = 0;
            accepted_steps = 0;
            count_total = 0;

        }
    }

    cout << energies*2*13.6 << endl;

    printFile(*file_energies, *file_energySquareds, *file_alpha, energies, energySquareds, alphas, betas);

}

double VMCSolver::localEnergy(const mat &r, const double &alpha, const double &beta)
{
    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    double waveFunctionCurrent = waveFunction(r, alpha, beta);

    // Kinetic energy

    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = waveFunction(rMinus, alpha, beta);
            waveFunctionPlus = waveFunction(rPlus, alpha, beta);
            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
    kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

    // Potential energy
    double potentialEnergy = 0;
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j)*r(i,j);
        }
        potentialEnergy -= charge / sqrt(rSingleParticle);
    }
    // Contribution from electron-electron potential
    double r12 = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            potentialEnergy += 1 / sqrt(r12);
        }
    }

    return kineticEnergy + potentialEnergy;
}



double VMCSolver::waveFunction(const mat &r, const double &alpha, const double &beta)
{
    double argument = 0;
    for(int i = 0; i < nParticles; i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        argument += sqrt(rSingleParticle);
    }

    double r12 = 0;
    double interaction = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            interaction += sqrt(r12)/(2*(1+beta*sqrt(r12)));
        }
    }


    return exp(-argument * alpha)*exp(interaction);
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
