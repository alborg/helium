#include "vmcsolver.h"

#include <iostream>
#include "../../../armadillo/include/armadillo"

using namespace arma;
using namespace std;

int main()
{
    double alpha = 0.0;
    double beta = 0.0;
    double alpha_step = 0.1;
    double beta_step = 0.1;
    int n_steps = 10;
    mat energies = zeros(n_steps,n_steps);
    mat energySquareds = zeros(n_steps,n_steps);


    VMCSolver *solver = new VMCSolver();

    for (int k=1; k<n_steps+1; k++) {
        alpha = k*alpha_step;
        for (int l=1; l<n_steps+1; l++) {
            beta = l*beta_step;
            mat energy = solver->runMonteCarloIntegration(alpha, beta);
            cout << alpha << " " << beta << energy << endl;
            //energies[k][l] = energy[0];
           // energySquareds[k][l] = energy[1];

        }
    }

    //cout << energies << endl;

    return 0;
}

