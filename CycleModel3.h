#ifndef CYCLEMODEL3_H
#define CYCLEMODEL3_H
#include "Matrice.h"
#include <string>

using namespace std;

Matrice CycleModel3(Matrice Pop, double mutT, double mutE, double mutA, double SigT, double SigE, double h, double S, double Rij, double Rjk, double Self, double I, double Concaveness, int Npop, string FIT, string ALLELES, string SHAPE);

#endif // CYCLEMODEL3_H
