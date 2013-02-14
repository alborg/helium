#include "vmcsolver.h"
#include "WaveFunction.h"
#include <armadillo>


using namespace arma;
using namespace std;

int main()
{
    //VMCSolver *solver = new VMCSolver();

    //solver->runMonteCarloIntegration();

    mat r("0.4 1.2 -0.3 ; -1.5 0.3 -0.5");
    cout << r << endl;
    WaveFunction *function = new WaveFunction(2, 3, 0.6, 0.9);
    double answer = function->waveFunction(r);
    cout << answer << endl;

    return 0;
}

