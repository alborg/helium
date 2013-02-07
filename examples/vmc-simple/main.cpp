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
    double alpha_min = 1.0;
    double alpha_max = 3.0;
    int alpha_steps = 2;
    double beta_min = 1.0;
    double beta_max = 3.0;
    int beta_steps = 2;

    double alpha_step = (alpha_max - alpha_min)/alpha_steps;
    double beta_step = (beta_max - beta_min)/beta_steps;

    mat energies = zeros(alpha_steps,beta_steps);
    mat energySquareds = zeros(alpha_steps,beta_steps);

    cout << alpha_steps << " " << beta_steps << endl;
    cout << alpha_step << " " << beta_step << endl;


    VMCSolver *solver = new VMCSolver();

    for (double k=alpha_min; k<alpha_max; k+=alpha_step) {
        alpha = k;
        for (double l=beta_min; l<beta_max; l+=beta_step) {
            beta = l;
            cout << alpha << " " << beta <<endl;
            mat energy = solver->runMonteCarloIntegration(alpha, beta);
            cout << "test" << endl;
            energies(k,l) = energy[0];
            energySquareds(k,l) = energy[1];


        }
    }

    cout << energies*2*13.6 << endl;

    ofstream myfile ("..\\..\\..\\output\\data.txt");
    if (myfile.is_open())
    {
        for (uint f=0; f<energies.n_rows; f++)
        {
            for (uint l=0; l<energies.n_cols; l++) {
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
        for (uint f=0; f<energies.n_rows; f++)
        {
            myfile2 << (f+1)*alpha_step << " ";
        }
        myfile2 << endl;
        for (uint l=0; l<energies.n_cols; l++) {
            myfile2 << (l+1)*beta_step << " ";
        }
        myfile2 << endl;


        myfile2.close();
    }
    else cout << "Unable to open file" << endl;


    return 0;
}

