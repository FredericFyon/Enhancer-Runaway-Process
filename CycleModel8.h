#ifndef CYCLEMODEL8_H_INCLUDED
#define CYCLEMODEL8_H_INCLUDED
#include "Matrice.h"
#include <string>

using namespace std;

Matrice CycleModel8(Matrice Pop, double Ei, double mutT, double mutE, double mutA, double SigT, double SigE, double h, double S, double I, double Rij, double Rjk, double Self, double Clone, int Npop, string FIT, string ALLELES, string SEX);

#endif // CYCLEMODEL8_H_INCLUDED
