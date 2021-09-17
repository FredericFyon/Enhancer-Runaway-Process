#ifndef CYCLEMODEL1_H
#define CYCLEMODEL1_H
#include "Matrice.h"
#include <string>

Matrice CycleModel1(Matrice Pop, double mutE, double mutA, double SigE, double h, double S, double Rij, double Rjk, double Self, double Clone, double Automixis, int Npop, string FIT, string ALLELES, string AUTOMIXIS);

#endif // CYCLEMODEL1_H
