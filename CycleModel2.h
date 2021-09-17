#ifndef CYCLEMODEL2_H
#define CYCLEMODEL2_H
#include <string>
#include "Matrice.h"

using namespace std;

Matrice CycleModel2(Matrice Pop, double mutE, double mutA, double SigE, double h, double S, double Rjk, double Self, double I, double Concaveness, int Npop, string FIT, string ALLELES, string SHAPE);

#endif // CYCLEMODEL2_H
