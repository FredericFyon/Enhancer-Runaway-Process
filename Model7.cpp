#include "Model7.h"
#include "CycleModel7.h"
#include "MoyFitModel7.h"
#include "Matrice.h"
#include "Statistics.h"
#include "LinearRegression.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>

using namespace std;

void Model7(double Ei, double mutE1, double mutE2, double mutA, double SigE1, double SigE2, double h, double S, double Rij, double Rjk, double Self, double I, int Npop, int Ngen, int Ngen2, int Nit, int TpsDesinit, int PasEchant, double Concaveness, string FIT, string ALLELES, string SUIVIFIT, string OutputPathway, string SHAPE)
{
    size_t TailleLineEdit = OutputPathway.size() + 1;
    char* FilePathway = new char[TailleLineEdit];
    strncpy(FilePathway,OutputPathway.c_str(),TailleLineEdit);                                  //Construction du fichier de sortie
    ofstream Results(FilePathway,ios::app);

    vector<double> E1Chrom1(Npop),E2Chrom1(Npop),SChrom1(Npop),E1Chrom2(Npop),E2Chrom2(Npop),SChrom2(Npop);
    vector<double> Empty;

    int indInitPop(0);
    for(indInitPop = 0 ; indInitPop < Npop ; ++indInitPop)
    {
        E1Chrom1[indInitPop] = Ei;
        E2Chrom1[indInitPop] = Ei;
        SChrom1[indInitPop] = 0.;
        E1Chrom2[indInitPop] = Ei;
        E2Chrom2[indInitPop] = Ei;
        SChrom2[indInitPop] = 0.;
    }
    Matrice PopInit(E1Chrom1,E2Chrom1,SChrom1,E1Chrom2,E2Chrom2,SChrom2,Empty,Empty,Empty,Empty);

    vector<vector<double> > VecLogE1It,VecLogE2It,VecWIt;

    double indIt(0);
    for(indIt = 0 ; indIt < Nit ; ++indIt)
    {
        //cout << "Iteration n°" << indIt+1 << " en cours" << endl;
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
            Pop = CycleModel7(Ei,Pop,0.,0.,mutA,SigE1,SigE2,h,S,Rij,Rjk,Self,I,Concaveness,Npop,FIT,ALLELES,SHAPE);
        }

        //cout << "Desinitialisation terminee" << endl;

        vector<double> VecLogE1,VecLogE2,VecW;

        int indGen(0);
        for(indGen = 0 ; indGen < Ngen ; ++indGen)
        {
            if(indGen%PasEchant == 0)
            {
                //cout << "Generation n°" << indGen << endl;
                vector<double> PopE1Chrom1 = Pop.getVec1();
                vector<double> PopE2Chrom1 = Pop.getVec2();
                vector<double> PopSChrom1 = Pop.getVec3();
                vector<double> PopE1Chrom2 = Pop.getVec4();
                vector<double> PopE2Chrom2 = Pop.getVec5();
                vector<double> PopSChrom2 = Pop.getVec6();

                vector<double> PopLogE1Chrom1,PopLogE2Chrom1,PopLogE1Chrom2,PopLogE2Chrom2;
                int indLog(0);
                for(indLog = 0 ; indLog < Npop ; ++indLog)
                {
                    PopLogE1Chrom1.push_back(log10(PopE1Chrom1[indLog]));
                    PopLogE2Chrom1.push_back(log10(PopE2Chrom1[indLog]));
                    PopLogE1Chrom2.push_back(log10(PopE1Chrom2[indLog]));
                    PopLogE2Chrom2.push_back(log10(PopE2Chrom2[indLog]));
                }

                Statistic PopLogE1(PopLogE1Chrom1,PopLogE1Chrom2);
                double PopLogE1Moy = PopLogE1.getMean();
                VecLogE1.push_back(PopLogE1Moy);

                Statistic PopLogE2(PopLogE2Chrom1,PopLogE2Chrom2);
                double PopLogE2Moy = PopLogE2.getMean();
                VecLogE2.push_back(PopLogE2Moy);

                if(SUIVIFIT == "Oui")
                {
                    double WMoy = MoyFitModel7(Ei,PopE1Chrom1,PopE2Chrom1,PopSChrom1,PopE1Chrom2,PopE2Chrom2,PopSChrom2,h,Npop,FIT,I);
                    VecW.push_back(WMoy);
                }
            }

            Pop = CycleModel7(Ei,Pop,0.,mutE2,mutA,SigE1,SigE2,h,S,Rij,Rjk,Self,I,Concaveness,Npop,FIT,ALLELES,SHAPE);
        }

        for(indGen = Ngen ; indGen < Ngen+Ngen2 ; ++indGen)
        {
            if(indGen%PasEchant == 0)
            {
                //cout << "Generation n°" << indGen << endl;
                vector<double> PopE1Chrom1 = Pop.getVec1();
                vector<double> PopE2Chrom1 = Pop.getVec2();
                vector<double> PopSChrom1 = Pop.getVec3();
                vector<double> PopE1Chrom2 = Pop.getVec4();
                vector<double> PopE2Chrom2 = Pop.getVec5();
                vector<double> PopSChrom2 = Pop.getVec6();

                vector<double> PopLogE1Chrom1,PopLogE2Chrom1,PopLogE1Chrom2,PopLogE2Chrom2;
                int indLog(0);
                for(indLog = 0 ; indLog < Npop ; ++indLog)
                {
                    PopLogE1Chrom1.push_back(log10(PopE1Chrom1[indLog]));
                    PopLogE2Chrom1.push_back(log10(PopE2Chrom1[indLog]));
                    PopLogE1Chrom2.push_back(log10(PopE1Chrom2[indLog]));
                    PopLogE2Chrom2.push_back(log10(PopE2Chrom2[indLog]));
                }

                Statistic PopLogE1(PopLogE1Chrom1,PopLogE1Chrom2);
                double PopLogE1Moy = PopLogE1.getMean();
                VecLogE1.push_back(PopLogE1Moy);

                Statistic PopLogE2(PopLogE2Chrom1,PopLogE2Chrom2);
                double PopLogE2Moy = PopLogE2.getMean();
                VecLogE2.push_back(PopLogE2Moy);

                if(SUIVIFIT == "Oui")
                {
                    double WMoy = MoyFitModel7(Ei,PopE1Chrom1,PopE2Chrom1,PopSChrom1,PopE1Chrom2,PopE2Chrom2,PopSChrom2,h,Npop,FIT,I);
                    VecW.push_back(WMoy);
                }
            }

            Pop = CycleModel7(Ei,Pop,mutE1,mutE2,mutA,SigE1,SigE2,h,S,Rij,Rjk,Self,I,Concaveness,Npop,FIT,ALLELES,SHAPE);
        }

        VecLogE1It.push_back(VecLogE1);
        VecLogE2It.push_back(VecLogE2);
        VecWIt.push_back(VecW);
    }

    //Results << "Fin simuls pop" << endl;

    /*int indVecX(0);
    vector<double> X;
    for(indVecX = 0 ; indVecX <= Ngen ; indVecX += PasEchant)
    {
        X.push_back(indVecX);                                                                   //Création du vecteur abscisse de la régression linéaire
    }

    int indPentes(0);
    vector<double> PentesE1,PentesE2;
    for(indPentes = 0 ; indPentes < Nit ; ++indPentes)
    {
        PentesE1.push_back(LinearRegression(X,VecLogE1It[indPentes]));
        PentesE2.push_back(LinearRegression(X,VecLogE2It[indPentes]));
    }*/

    vector<double> EvolLogE1((Ngen+Ngen2)/PasEchant+1),EvolLogE2((Ngen+Ngen2)/PasEchant+1),EvolW((Ngen+Ngen2)/PasEchant+1);

    if(SUIVIFIT == "Oui")
    {
        int indMoy(0);
        for(indMoy = 0 ; indMoy < (Ngen+Ngen2)/PasEchant ; ++indMoy)
        {
            EvolLogE1[indMoy] = 0;
            EvolLogE2[indMoy] = 0;
            EvolW[indMoy] = 0;
            int indTampon(0);
            for(indTampon ; indTampon < Nit ; ++indTampon)
            {
                EvolLogE1[indMoy] = EvolLogE1[indMoy] + VecLogE1It[indTampon][indMoy];
                EvolLogE2[indMoy] = EvolLogE2[indMoy] + VecLogE2It[indTampon][indMoy];
                EvolW[indMoy] = EvolW[indMoy] + VecWIt[indTampon][indMoy];
            }
            EvolLogE1[indMoy] = EvolLogE1[indMoy]/Nit;
            EvolLogE2[indMoy] = EvolLogE2[indMoy]/Nit;
            EvolW[indMoy] = EvolW[indMoy]/Nit;
        }
    }
    else if(SUIVIFIT == "Non")
    {
        int indMoy(0);
        for(indMoy = 0 ; indMoy < (Ngen+Ngen2)/PasEchant ; ++indMoy)
        {
            EvolLogE1[indMoy] = 0;
            EvolLogE2[indMoy] = 0;
            int indTampon(0);
            for(indTampon ; indTampon < Nit ; ++indTampon)
            {
                EvolLogE1[indMoy] = EvolLogE1[indMoy] + VecLogE1It[indTampon][indMoy];
                EvolLogE2[indMoy] = EvolLogE2[indMoy] + VecLogE2It[indTampon][indMoy];
            }
            EvolLogE1[indMoy] = EvolLogE1[indMoy]/Nit;
            EvolLogE2[indMoy] = EvolLogE2[indMoy]/Nit;
        }
    }

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
        Results << "Taux de mutation au Locus E2 : " << mutE2 << endl;
        Results << "Taux de mutation au Locus E1 : " << mutE1 << endl;
        Results << "Ecart-type de la loi normale régissant les mutations multiplicatives au Locus E2 : " << SigE2 << endl;
        Results << "Ecart-type de la loi normale régissant les mutations multiplicatives au Locus E1 : " << SigE1 << endl;
        Results << "Taux de recombinaison entre les locus T et E : " << Rij << endl;
        Results << "Taux de recombinaison entre les locus E et A : " << Rjk << endl;
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
        /*Results << "Pentes des régressions linéaires du log des forces des Cis-acting factors en fonction du temps : " << endl;
        Results << "LocusE1     LocusE2" << endl;
        int indPente(0);
        for(indPente = 0 ; indPente < Nit ; ++indPente)
        {
            Results.precision(10);
            Results << fixed;
            Results << PentesE1[indPente] << " " << PentesE2[indPente] << endl;
        }
        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;*/
        Results << "Evolution Moyenne" << endl;
        Results << endl;

        if(SUIVIFIT == "Non")
        {
            int indRes(0);
            for(indRes = 0 ; indRes < (Ngen+Ngen2)/PasEchant; ++indRes)
            {
                Results << indRes*PasEchant << " " << EvolLogE1[indRes] << " " << EvolLogE2[indRes] << endl;
            }
            Results << endl;
        }

        else if(SUIVIFIT == "Oui")
        {
            int indRes(0);
            for(indRes = 0 ; indRes < (Ngen+Ngen2)/PasEchant ; ++indRes)
            {
                Results << indRes*PasEchant << " " << EvolLogE1[indRes] << " " << EvolLogE2[indRes] << " " << EvolW[indRes] << endl;
            }
            Results << endl;
        }

        if(SUIVIFIT == "Non")
        {
            Results << "Moyennes dans les populations du log des forces des Cis-acting factors : " << endl;
            Results << "Générations    LocusE1    LocusE2" << endl;

            int indIteration(0);
            for(indIteration = 0 ; indIteration < Nit ; ++ indIteration)
            {
                Results << endl;
                Results << "Itération n°" << indIteration+1 << endl;

                int indGeneration(0);
                for(indGeneration = 0 ; indGeneration <= (Ngen+Ngen2)/PasEchant ; ++indGeneration)
                {
                    Results.precision(6);
                    Results << fixed;
                    Results << indGeneration*PasEchant << " " << VecLogE1It[indIteration][indGeneration] << " " << VecLogE2It[indIteration][indGeneration] << endl;
                }
            }
        }
        else if(SUIVIFIT == "Oui")
        {
            Results << "Moyennes dans les populations du log des forces des Cis-acting factors et des fitness : " << endl;
            Results << "Générations    LocusE1    LocusE2   Fitness" << endl;

            int indIteration(0);
            for(indIteration = 0 ; indIteration < Nit ; ++ indIteration)
            {
                Results << endl;
                Results << "Itération n°" << indIteration+1 << endl;

                int indGeneration(0);
                for(indGeneration = 0 ; indGeneration <= (Ngen+Ngen2)/PasEchant ; ++indGeneration)
                {
                    Results.precision(6);
                    Results << fixed;
                    Results << indGeneration*PasEchant << " " << VecLogE1It[indIteration][indGeneration] << " " << VecLogE2It[indIteration][indGeneration] << " " << VecWIt[indIteration][indGeneration] << endl;
                }
            }
        }
    }
}
