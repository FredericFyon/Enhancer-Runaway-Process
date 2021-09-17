#include "LinearRegression.h"
#include "Statistics.h"
#include <math.h>
#include <vector>

using namespace std;

double LinearRegression(vector<double> x, vector<double> y)
{
    double Taille = x.size();

    Statistic X(x,x);
    double xMoy = X.getMean();
    double xVar = X.getVariance();

    Statistic Y(y,y);
    double yMoy = Y.getMean();

    double cov(0);
    int i(0);

    for(i = 0 ; i < Taille ; i++)
    {
        cov += (x[i] - xMoy)*(y[i] - yMoy);
    }

    cov /= Taille;

    double pente = cov/xVar;

    return pente;
}
