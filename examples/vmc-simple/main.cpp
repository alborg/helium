#include "vmcsolver.h"

#include <armadillo>


using namespace arma;
using namespace std;

int main()
{
    VMCSolver *solver = new VMCSolver();

    solver->runMonteCarloIntegration();

    return 0;
}

