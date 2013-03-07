#include "vmcsolver.h"
#include "vmcimportance.h"
#include <armadillo>


using namespace arma;
using namespace std;

int main(int argc, char* argv[])
{

    //VMCSolver *solver = new VMCSolver();
    VMCImportance *importanceSolver = new VMCImportance();

    //solver->runMonteCarloIntegration(argc, argv);
    importanceSolver->runMonteCarloIntegration(argc, argv);


    return 0;
}

