#ifndef CYCLEMODEL6_H_INCLUDED
#define CYCLEMODEL6_H_INCLUDED
#include <string>
#include "Matrice.h"

using namespace std;

Matrice CycleModel6(Matrice Pop, double Ei, double mutD, double mutE, double mutA, double SigD, double SigE, double h, double S, double Rij, double Rjk, double Self, double I, int Npop, string FIT, string ALLELES, string Modifier);

#endif // CYCLEMODEL6_H_INCLUDED
