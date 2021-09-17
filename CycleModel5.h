#ifndef CYCLEMODEL5_H_INCLUDED
#define CYCLEMODEL5_H_INCLUDED
#include "Matrice.h"
#include <string>

using namespace std;

Matrice CycleModel5(Matrice Pop, double Ei, double mutD, double mutE, double mutA, double SigD, double SigE, double h, double S, double Rij, double Rjk, double Self, double I, int Npop, string FIT, string ALLELES);


#endif // CYCLEMODEL5_H_INCLUDED
