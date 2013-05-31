#include "vmcsolver.h"
#include <armadillo>


using namespace arma;
using namespace std;

int main(int argc, char* argv[])
{

    //Start the solver:
    VMCSolver *solver = new VMCSolver();
    solver->runMonteCarloIntegration(argc, argv);


    return 0;
}

