#include "CalculateVariance.h"
#include <vector>

using namespace std;

double CalculateVariance(vector<double> Vec, double mean)
{
    int n = Vec.size();
    double Var(0);
    int i = 0;

    for(i = 0 ; i < n ; ++i)
    {
        Var = Var + (Vec[i] - mean)*(Vec[i] - mean);
    }

    Var = Var / n;

    return Var;
}
