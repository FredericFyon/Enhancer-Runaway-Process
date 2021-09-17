#ifndef CYCLEMODEL7_H_INCLUDED
#define CYCLEMODEL7_H_INCLUDED

#include <string>
#include "Matrice.h"

using namespace std;

Matrice CycleModel7(double Ei, Matrice Pop, double mutE1, double mutE2, double mutA, double SigE1, double SigE2, double h, double S, double Rij, double Rjk, double Self, double I, double Concaveness, int Npop, string FIT, string ALLELES, string SHAPE);


#endif // CYCLEMODEL7_H_INCLUDED
