#ifndef CYCLEMODEL4_H_INCLUDED
#define CYCLEMODEL4_H_INCLUDED
#include <string>
#include "Matrice.h"

using namespace std;

Matrice CycleModel4(Matrice Pop, double Ei, double mutR, double mutE, double mutA, double SigR, double SigE, double h, double S, double Rij, double Self, double I, int Npop, string FIT, string ALLELES);


#endif // CYCLEMODEL4_H_INCLUDED
