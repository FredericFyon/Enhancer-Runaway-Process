#include "CalculateMean.h"
#include <vector>

using namespace std;

double CalculateMean(vector<double> Vec)
{
    int n = Vec.size();
    double mean(0);
    int i(0);
    for(i = 0 ; i < n ; ++i)
    {
        mean = mean + Vec[i];
    }
    mean = mean/n;

    return mean;
}
