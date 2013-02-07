#include "vmcsolver.h"

#include <iostream>
#include <fstream>
#include <armadillo>


using namespace arma;
using namespace std;

int main()
{
    double alpha = 0.0;
    double beta = 0.0;
    double alpha_step = 0.5;
    double beta_step = 0.5;
    int n_steps = 5;
    mat energies = zeros(n_steps,n_steps);
    mat energySquareds = zeros(n_steps,n_steps);


    VMCSolver *solver = new VMCSolver();

    for (int k=0; k<n_steps; k++) {
        alpha = (k+1)*alpha_step;
        for (int l=0; l<n_steps; l++) {
            beta = (l+1)*beta_step;
            mat energy = solver->runMonteCarloIntegration(alpha, beta);
            energies(k,l) = energy[0];
            energySquareds(k,l) = energy[1];

        }
    }

    cout << energies << endl;

    ofstream myfile ("..\\..\\..\\output\\data.txt");
    if (myfile.is_open())
    {
        for (int f=0; f<energies.n_rows; f++)
        {
            for (int l=0; l<energies.n_cols; l++) {
                myfile << energies(f,l)*2*13.6 << " ";
            }
            myfile << endl;
        }

        myfile.close();
    }
    else cout << "Unable to open file" << endl;


    ofstream myfile2 ("..\\..\\..\\output\\data_alpha_beta.txt");
    if (myfile2.is_open())
    {
        for (int f=0; f<energies.n_rows; f++)
        {
            myfile2 << (f+1)*alpha_step << " ";
        }
        myfile2 << endl;
        for (int l=0; l<energies.n_cols; l++) {
            myfile2 << (l+1)*beta_step << " ";
        }
        myfile2 << endl;


        myfile2.close();
    }
    else cout << "Unable to open file" << endl;


    return 0;
}

