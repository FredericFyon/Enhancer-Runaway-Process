#include "Model2.h"
#include "CycleModel2.h"
#include "MoyFitModel2.h"
#include "Matrice.h"
#include "Statistics.h"
#include "LinearRegression.h"
#include <fstream>
#include <cmath>
#include <string>
#include <cstring>

using namespace std;

void Model2(double Ei, double mutE, double mutA, double SigE, double h, double S, double Rjk, double Self, double Concaveness, double I, int Npop, int Ngen, int Nit, int TpsDesinit, int PasEchant, string FIT, string ALLELES, string SUIVIFIT, string OutputPathway, string SHAPE)
{
    size_t TailleLineEdit = OutputPathway.size() + 1;
    char* FilePathway = new char[TailleLineEdit];
    strncpy(FilePathway,OutputPathway.c_str(),TailleLineEdit);                                  //Construction du fichier de sortie
    ofstream Results(FilePathway,ios::app);

    vector<double> EChrom1a(Npop), SChrom1a(Npop), EChrom2a(Npop), SChrom2a(Npop), EChrom1b(Npop), SChrom1b(Npop), EChrom2b(Npop), SChrom2b(Npop);
    vector<double> Empty;

    int indInitPop(0);
    for(indInitPop = 0 ; indInitPop < Npop ; ++indInitPop)
    {
        EChrom1a[indInitPop] = Ei;
        SChrom1a[indInitPop] = 0.;
        EChrom2a[indInitPop] = Ei;
        SChrom2a[indInitPop] = 0.;
        EChrom1b[indInitPop] = Ei;
        SChrom1b[indInitPop] = 0.;
        EChrom2b[indInitPop] = Ei;
        SChrom2b[indInitPop] = 0.;
    }
    Matrice PopInit(EChrom1a,SChrom1a,EChrom2a,SChrom2a,EChrom1b,SChrom1b,EChrom2b,SChrom2b,Empty,Empty);

    vector<vector<double> > VecLogEaIt,VecLogEbIt,VecWIt;

    double indIt(0);
    for(indIt = 0 ; indIt < Nit ; ++indIt)
    {
        if(Results)
        {
            Results << "Iteration n°" << indIt+1 << " in progress..." << endl;
        }

        Matrice Pop(Empty,Empty,Empty,Empty,Empty,Empty,Empty,Empty,Empty,Empty);

        int indDesinit(0);
        for(indDesinit = 0 ; indDesinit < TpsDesinit ; ++indDesinit)
        {
            if(indDesinit == 0)
            {
                Pop = PopInit;
            }
            Pop = CycleModel2(Pop,0.,mutA,SigE,h,S,Rjk,Self,I,Concaveness,Npop,FIT,ALLELES,SHAPE);
        }

        vector<double> VecLogEa,VecLogEb,VecW;

        int indGen(0);
        for(indGen = 0 ; indGen <= Ngen ; ++indGen)
        {
            if(indGen%PasEchant == 0)
            {
                vector<double> PopEChrom1a = Pop.getVec1();
                vector<double> PopSChrom1a = Pop.getVec2();
                vector<double> PopEChrom2a = Pop.getVec3();
                vector<double> PopSChrom2a = Pop.getVec4();
                vector<double> PopEChrom1b = Pop.getVec5();
                vector<double> PopSChrom1b = Pop.getVec6();
                vector<double> PopEChrom2b = Pop.getVec7();
                vector<double> PopSChrom2b = Pop.getVec8();

                vector<double> PopLogEChrom1a,PopLogEChrom2a,PopLogEChrom1b,PopLogEChrom2b;
                int indLog(0);
                for(indLog = 0 ; indLog < PopEChrom1a.size() ; ++indLog)
                {
                    PopLogEChrom1a.push_back(log10(PopEChrom1a[indLog]));
                    PopLogEChrom2a.push_back(log10(PopEChrom2a[indLog]));
                    PopLogEChrom1b.push_back(log10(PopEChrom1b[indLog]));
                    PopLogEChrom2b.push_back(log10(PopEChrom2b[indLog]));
                }

                Statistic PopLogEa(PopLogEChrom1a,PopLogEChrom2a);
                double PopLogEaMoy = PopLogEa.getMean();
                VecLogEa.push_back(PopLogEaMoy);

                Statistic PopLogEb(PopLogEChrom1b,PopLogEChrom2b);
                double PopLogEbMoy = PopLogEb.getMean();
                VecLogEb.push_back(PopLogEbMoy);

                if(SUIVIFIT == "Oui")
                {
                    double WMoy = MoyFitModel2(PopEChrom1a,PopSChrom1a,PopEChrom2a,PopSChrom2a,PopEChrom1b,PopSChrom1b,PopEChrom2b,PopSChrom2b,h,Npop,I,Concaveness,FIT,SHAPE);
                    VecW.push_back(WMoy);
                }
            }

            Pop = CycleModel2(Pop,mutE,mutA,SigE,h,S,Rjk,Self,I,Concaveness,Npop,FIT,ALLELES,SHAPE);
        }

        VecLogEaIt.push_back(VecLogEa);
        VecLogEbIt.push_back(VecLogEb);
        VecWIt.push_back(VecW);
    }

    vector<double> VecMoyLogEa(Ngen/PasEchant+1), VecMoyLogEb(Ngen/PasEchant+1);
    int indMoy(0);
    for(indMoy = 0 ; indMoy <= Ngen/PasEchant ; ++indMoy)
    {
        VecMoyLogEa[indMoy] = 0;
        VecMoyLogEb[indMoy] = 0;
        int indTaMere(0);
        for(indTaMere = 0 ; indTaMere < Nit ; ++indTaMere)
        {
            VecMoyLogEa[indMoy] = VecMoyLogEa[indMoy] + VecLogEaIt[indTaMere][indMoy];
            VecMoyLogEb[indMoy] = VecMoyLogEb[indMoy] + VecLogEbIt[indTaMere][indMoy];
        }
        VecMoyLogEa[indMoy] = VecMoyLogEa[indMoy]/Nit;
        VecMoyLogEb[indMoy] = VecMoyLogEb[indMoy]/Nit;
    }

    int indVecX(0);
    vector<double> X;
    for(indVecX = 0 ; indVecX <= Ngen ; indVecX += PasEchant)
    {
        X.push_back(indVecX);                                                                   //Création du vecteur abscisse de la régression linéaire
    }

    int indPentes(0);
    double PentesA,PentesB;
    PentesA = LinearRegression(X,VecMoyLogEa);
    PentesB = LinearRegression(X,VecMoyLogEb);


    if(Results)
    {
        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        Results << "RESULTATS" << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        Results << "Valeurs des paramètres : " << endl;
        Results << endl;
        Results << "Force initiale des Cis-acting factors : " << Ei << endl;
        Results << "Coefficient de dominance : " << h << endl;
        Results << "Intensité de la sélection au Locus A : " << S << endl;
        Results << "Taux de mutation au Locus A : " << mutA << endl;
        Results << "Taux de mutation au Locus E : " << mutE << endl;
        Results << "Ecart-type de la loi normale régissant les mutations multiplicatives au Locus E : " << SigE << endl;
        Results << "Taux de recombinaison entre les deux locus : " << Rjk << endl;
        Results << "Taux d'autofécondation : " << Self << endl;
        Results << "Intensité de la contrainte évolutive : " << I << endl;
        Results << "Concavité de la contrainte évolutive : " << Concaveness << endl;
        Results << "Taille de la population d'individus diploïdes : " << Npop << endl;
        Results << "Nombre de générations simulées : " << Ngen << endl;
        Results << "Nombre d'itérations du processus évolutif simulées : " << Nit << endl;
        Results << "Nombre de générations simulées pour désinitialiser la population : " << TpsDesinit << endl;
        Results << "Pas d'échantillonage des valeurs génomiques : " << PasEchant << endl;
        if(FIT == "Selection")
        {
            Results << "Sélection active" << endl;
        }
        else if(FIT == "Derive")
        {
            Results << "Dérive" << endl;
        }
        if(ALLELES == "2")
        {
            Results << "2 allèles au Locus A" << endl;
        }
        else if(ALLELES == "Infinite")
        {
            Results << "Nombre infini d'allèles au locus A" << endl;
        }
        if(SHAPE == "Linear")
        {
            Results << "Contrainte évolutive linéaire" << endl;
        }
        else if(SHAPE == "Concave")
        {
            Results << "Contrainte évolutive concave" << endl;
        }
        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        Results << "Pentes des régressions linéaires du log des forces des Cis-acting factors en fonction du temps : " << endl;
        Results << "1ère paire     2ème paire" << endl;
        Results.precision(10);
        Results << fixed;
        Results << PentesA << " " << PentesB << endl;
        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        if(SUIVIFIT == "Non")
        {
            Results << "Moyennes dans les populations du log des forces des Cis-acting factors : " << endl;
            Results << "Générations    1ère paire    2ème paire" << endl;

            int indIteration(0);
            for(indIteration = 0 ; indIteration < Nit ; ++ indIteration)
            {
                Results << endl;
                Results << "Itération n°" << indIteration+1 << endl;

                int indGeneration(0);
                for(indGeneration = 0 ; indGeneration <= Ngen/PasEchant ; ++indGeneration)
                {
                    Results.precision(6);
                    Results << fixed;
                    Results << indGeneration*PasEchant << " " << VecLogEaIt[indIteration][indGeneration] << " " << VecLogEbIt[indIteration][indGeneration] << endl;
                }
            }
        }
        else if(SUIVIFIT == "Oui")
        {
            Results << "Moyennes dans les populations du log des forces des Cis-acting factors et des fitness : " << endl;
            Results << "Générations    1ère paire    2ème paire   Fitness" << endl;

            int indIteration(0);
            for(indIteration = 0 ; indIteration < Nit ; ++ indIteration)
            {
                Results << endl;
                Results << "Itération n°" << indIteration+1 << endl;

                int indGeneration(0);
                for(indGeneration = 0 ; indGeneration <= Ngen/PasEchant ; ++indGeneration)
                {
                    Results.precision(6);
                    Results << fixed;
                    Results << indGeneration*PasEchant << " " << VecLogEaIt[indIteration][indGeneration] << " " << VecLogEbIt[indIteration][indGeneration] << " " << VecWIt[indIteration][indGeneration] << endl;
                }
            }
        }
    }
}

