#include "Model8.h"
#include "CycleModel8.h"
#include "MoyFitModel1.h"
#include "Matrice.h"
#include "Statistics.h"
#include "LinearRegression.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <cstring>
#include "CalculateMean.h"
#include "CalculateVariance.h"
#include <boost/random.hpp>


using namespace std;

void Model8(double Ei, double mutT, double mutE, double mutA, double SigT, double SigE, double h, double S, double I, double Rij, double Rjk, double Self, double Clone, int Npop, int Ngen, int Nit, int TpsDesinit, int PasEchant, string FIT, string ALLELES, string SUIVIFIT, string OutputPathway, string Log, string SEX)
{
    size_t TailleLineEdit = OutputPathway.size() + 1;
    char* FilePathway = new char[TailleLineEdit];
    strncpy(FilePathway,OutputPathway.c_str(),TailleLineEdit);                                  //Construction du fichier de sortie
    ofstream Results(FilePathway,ios::app);

    vector<double> SexChrom1(Npop),TChrom1(Npop),EChrom1(Npop),SChrom1(Npop),SexChrom2(Npop),TChrom2(Npop),EChrom2(Npop),SChrom2(Npop);                     //Création des 4 vecteurs représentant les 4 valeurs génomiques de chaque individu
    vector<double> Empty;                                                                       //Création d'un vecteur vide

    int indInitPop(0);
    for(indInitPop = 0 ; indInitPop < Npop ; ++indInitPop)
    {
        EChrom1[indInitPop] = Ei;
        EChrom2[indInitPop] = Ei;
        TChrom1[indInitPop] = 10.;
        TChrom2[indInitPop] = 10.;
        SChrom1[indInitPop] = 0.;                                                               //Initiation des valeurs génomiques : population monomorphique
        SChrom2[indInitPop] = 0.;
        if(indInitPop < Npop/2)
        {
            SexChrom1[indInitPop] = 0;
            SexChrom2[indInitPop] = 0;
        }
        else
        {
            double chrom = rand() %2 + 1;
            if(chrom == 1)
            {
                SexChrom1[indInitPop] = 0;
                SexChrom2[indInitPop] = 1;
            }
            else if(chrom == 2)
            {
                SexChrom1[indInitPop] = 1;
                SexChrom2[indInitPop] = 0;
            }
        }
    }
    Matrice PopInit(SexChrom1,TChrom1,EChrom1,SChrom1,SexChrom2,TChrom2,EChrom2,SChrom2,Empty,Empty);       //Matrice représentant la population initiale

    vector<vector<double> > VecLogTIt,VecLogEIt,VecWIt,VarASEVecIt,VecEIt,VecHIt,VecLog0It,VecLog1It,VecLogS1It,VecLogS0It,VecS0It,VecS1It,VecXmaleIt,VecYIt,VecXfem1It,VecXfem2It;
    vector<vector<vector<double> > > VecASEIt;                                                   //Création des matrices stockant les forces des Cis-acting factors et les fitness au cours
                                                                                                //des générations et pour les Nit itérations
    double indIt(0);
    for(indIt = 0 ; indIt < Nit ; ++indIt)                                                      //Boucle des Nit itérations de l'évolution d'une même population initiale
    {
        vector<double> VecLogT,VecLogE,VecW,VarASEVec,VecE,VecH,VecLog0,VecLog1,VecLogS0,VecLogS1,VecS0,VecS1,VecXmale,VecY,VecXfem1,VecXfem2;
        vector<vector<double> > VecASETps;

        if(Results)
        {
            Results << "Iteration n°" << indIt+1 << " in progress..." << endl;
        }

        Matrice Pop(Empty,Empty,Empty,Empty,Empty,Empty,Empty,Empty,Empty,Empty);               //Création d'une matrice représentant la population

        int indDesinit(0);
        for(indDesinit = 0 ; indDesinit < TpsDesinit ; ++indDesinit)                            //TpsDesinit générations sont simulées pour atteindre l'équilibre mutation-sélection
        {                                                                                       //pour la fréquence en allèle délétère du gène. Pas de mutation au locus du Cis-acting factor
            //cout << indDesinit << endl;
            if(indDesinit == 0)
            {
                Pop = PopInit;                                                                  //Avant que la désinitialisation commence, la population est la population initiale
            }
            Pop = CycleModel8(Pop,Ei,0.,0.,mutA,SigT,SigE,h,S,I,Rij,Rjk,Self,Clone,Npop,FIT,ALLELES,SEX);
            /*if(indDesinit%PasEchant == 0)
            {
                vector<double> PopH = Pop.getVec7();
                Statistic PopHMOY(PopH,PopH);
                double PopHMoy = PopHMOY.getMean();
                VecH.push_back(PopHMoy);
                if(SUIVIFIT == "Oui")
                {
                    double WMoy = MoyFitModel1(Pop.getVec2(),Pop.getVec3(),Pop.getVec5(),Pop.getVec6(),h,Npop,FIT);    //Stockage à intervalles réguliers de la moyenne de la fitness dans la population
                    //Results << WMoy << endl;
                    VecW.push_back(WMoy);
                }
            }*/
        }

        //cout << "Fin Desinitialisation" << endl;
                                                                //Création des vecteurs stockant les forces des Cis-acting factors et les fitness au cours
                                                                                                //des générations
        int indGen(0);
        for(indGen = 0 ; indGen <= Ngen ; ++indGen)
        {
            //cout << indIt << "; " << indGen << endl;
            if(indGen%PasEchant == 0)
            {
                //cout << indGen << endl;
                vector<double> PopSexChrom1 = Pop.getVec1();
                vector<double> PopTChrom1 = Pop.getVec2();
                vector<double> PopEChrom1 = Pop.getVec3();
                vector<double> PopSChrom1 = Pop.getVec4();
                vector<double> PopSexChrom2 = Pop.getVec5();
                vector<double> PopTChrom2 = Pop.getVec6();
                vector<double> PopEChrom2 = Pop.getVec7();                                      //Extraction à intervalle régulier des valeurs génomiques de la population
                vector<double> PopSChrom2 = Pop.getVec8();
                vector<double> PopH = Pop.getVec9();
                vector<double> PopExpr = Pop.getVec10();
                //cout << PopExpr[0] << " " << PopExpr[1] << " " << PopExpr[2] << " " << PopExpr[3] << endl;

                vector<double> PopLogTChrom1,PopLogEChrom1,PopLogTChrom2,PopLogEChrom2,PopLog0,PopLog1,PopLogS0,PopLogS1,PopS0,PopS1;
                int indLog(0);
                for(indLog = 0 ; indLog < Npop ; ++indLog)
                {
                    PopLogTChrom1.push_back(log10(PopTChrom1[indLog]));
                    PopLogTChrom2.push_back(log10(PopTChrom2[indLog]));
                    PopLogEChrom1.push_back(log10(PopEChrom1[indLog]));                         //Passage des données au logarithme décimal
                    PopLogEChrom2.push_back(log10(PopEChrom2[indLog]));
                    if(PopSexChrom1[indLog] == 0)
                    {
                        PopLog0.push_back(log10(PopEChrom1[indLog]));
                        PopS0.push_back(PopSChrom1[indLog]);
                        PopLogS0.push_back(log10(PopSChrom1[indLog]));
                    }
                    else if(PopSexChrom1[indLog] == 1)
                    {
                        PopLog1.push_back(log10(PopEChrom1[indLog]));
                        PopS1.push_back(PopSChrom1[indLog]);
                        PopLogS1.push_back(log10(PopSChrom1[indLog]));
                    }

                    if(PopSexChrom2[indLog] == 0)
                    {
                        PopLog0.push_back(log10(PopEChrom2[indLog]));
                        PopS0.push_back(PopSChrom2[indLog]);
                        PopLogS0.push_back(log10(PopSChrom2[indLog]));
                    }
                    else if(PopSexChrom2[indLog] == 1)
                    {
                        PopLog1.push_back(log10(PopEChrom2[indLog]));
                        PopS1.push_back(PopSChrom2[indLog]);
                        PopLogS1.push_back(log10(PopSChrom2[indLog]));
                    }
                }

                Statistic PopLogE(PopLogEChrom1,PopLogEChrom2);
                double PopLogEMoy = PopLogE.getMean();
                VecLogE.push_back(PopLogEMoy);
                Statistic PopLogT(PopLogTChrom1,PopLogTChrom2);
                double PopLogTMoy = PopLogT.getMean();
                VecLogT.push_back(PopLogTMoy);
                Statistic PopE(PopEChrom1,PopEChrom2);
                double PopEMoy = log10(PopE.getMean());
                VecE.push_back(PopEMoy);
                Statistic PopHMOY(PopH,PopH);
                double PopHMoy = PopHMOY.getMean();
                VecH.push_back(PopHMoy);
                Statistic PopLogSex0(PopLog0,PopLog0);
                double PopLog0Moy = PopLogSex0.getMean();
                VecLog0.push_back(PopLog0Moy);
                Statistic PopLogSex1(PopLog1,PopLog1);
                double PopLog1Moy = PopLogSex1.getMean();
                VecLog1.push_back(PopLog1Moy);
                Statistic PopLogSexS0(PopLogS0,PopLogS0);
                double PopLogS0Moy = PopLogSexS0.getMean();
                VecLogS0.push_back(PopLogS0Moy);
                Statistic PopLogSexS1(PopLogS1,PopLogS1);
                double PopLogS1Moy = PopLogSexS1.getMean();
                //cout << PopLogS1Moy << endl;
                Statistic PopSexS0(PopS0,PopS0);
                double PopS0Moy = PopSexS0.getMean();
                VecS0.push_back(PopS0Moy);
                Statistic PopSexS1(PopS1,PopS1);
                double PopS1Moy = PopSexS1.getMean();
                VecLogS1.push_back(PopLogS1Moy);
                VecS1.push_back(PopS1Moy);
                VecXmale.push_back(PopExpr[0]);
                VecY.push_back(PopExpr[1]);
                VecXfem1.push_back(PopExpr[2]);
                VecXfem2.push_back(PopExpr[3]);
                                                                //Stockage à intervalles réguliers de la moyenne du log des forces des Cis-acting factors
                if(SUIVIFIT == "Oui")
                {
                    double WMoy = MoyFitModel1(PopEChrom1,PopSChrom1,PopEChrom2,PopSChrom2,h,Npop,FIT);    //Stockage à intervalles réguliers de la moyenne de la fitness dans la population
                    //Results << WMoy << endl;
                    VecW.push_back(WMoy);
                }
            }

            Pop = CycleModel8(Pop,Ei,mutT,mutE,mutA,SigT,SigE,h,S,I,Rij,Rjk,Self,Clone,Npop,FIT,ALLELES,SEX);
        }

        VecLogEIt.push_back(VecLogE);
        VecLogTIt.push_back(VecLogT);
        VecLog0It.push_back(VecLog0);
        VecLog1It.push_back(VecLog1);
        VecLogS0It.push_back(VecLogS0);
        VecLogS1It.push_back(VecLogS1);
        VecS0It.push_back(VecS0);
        VecS1It.push_back(VecS1);
        VecEIt.push_back(VecE);
        VecWIt.push_back(VecW);
        VecHIt.push_back(VecH);                                                               //Stockage à chaque itération de l'évolution de la moyenne dans la pop du log des forces des
        VecASEIt.push_back(VecASETps);
        VarASEVecIt.push_back(VarASEVec);
        VecXmaleIt.push_back(VecXmale);
        VecYIt.push_back(VecY);
        VecXfem1It.push_back(VecXfem1);
        VecXfem2It.push_back(VecXfem2);                                                                                     //Cis-acting factors et de la moyenne dans la pop des fitness
    }

    vector<double> VecMoyXmale(Ngen/PasEchant+1),VecMoyY(Ngen/PasEchant+1),VecMoyXfem1(Ngen/PasEchant+1),VecMoyXfem2(Ngen/PasEchant+1),VecMoyLogT(Ngen/PasEchant+1),VecMoyLogE(Ngen/PasEchant+1), VecMoyE(Ngen/PasEchant+1),VecMoyLog0(Ngen/PasEchant+1),VecMoyLog1(Ngen/PasEchant+1),VecMoyS0(Ngen/PasEchant+1),VecMoyS1(Ngen/PasEchant+1),VecMoyH((Ngen)/PasEchant+1),VecMoyW((Ngen+TpsDesinit)/PasEchant + 1);
    int indMoy(0);
    for(indMoy = 0 ; indMoy <= Ngen/PasEchant ; ++indMoy)
    {
        VecMoyLogE[indMoy] = 0;
        VecMoyLogT[indMoy] = 0;
        VecMoyE[indMoy] = 0;
        VecMoyLog0[indMoy] = 0;
        VecMoyLog1[indMoy] = 0;
        VecMoyS0[indMoy] = 0;
        VecMoyS1[indMoy] = 0;
        VecMoyXmale[indMoy] = 0;
        VecMoyY[indMoy] = 0;
        VecMoyXfem1[indMoy] = 0;
        VecMoyXfem2[indMoy] = 0;
        VecMoyH[indMoy] = 0;

        int indTaMere(0);
        for(indTaMere = 0 ; indTaMere < Nit ; ++indTaMere)
        {
            VecMoyLogE[indMoy] = VecMoyLogE[indMoy] + VecLogEIt[indTaMere][indMoy];
            VecMoyLogT[indMoy] = VecMoyLogT[indMoy] + VecLogTIt[indTaMere][indMoy];
            VecMoyE[indMoy] = VecMoyE[indMoy] + VecEIt[indTaMere][indMoy];
            VecMoyLog0[indMoy] = VecMoyLog0[indMoy] + VecLog0It[indTaMere][indMoy];
            VecMoyLog1[indMoy] = VecMoyLog1[indMoy] + VecLog1It[indTaMere][indMoy];
            VecMoyS0[indMoy] = VecMoyS0[indMoy] + VecS0It[indTaMere][indMoy];
            VecMoyS1[indMoy] = VecMoyS1[indMoy] + VecS1It[indTaMere][indMoy];
            VecMoyXmale[indMoy] = VecMoyXmale[indMoy] + VecXmaleIt[indTaMere][indMoy];
            VecMoyY[indMoy] = VecMoyY[indMoy] + VecYIt[indTaMere][indMoy];
            VecMoyXfem1[indMoy] = VecMoyXfem1[indMoy] + VecXfem1It[indTaMere][indMoy];
            VecMoyXfem2[indMoy] = VecMoyXfem2[indMoy] + VecXfem2It[indTaMere][indMoy];
            VecMoyH[indMoy] = VecMoyH[indMoy] + VecHIt[indTaMere][indMoy];
        }
        VecMoyLogE[indMoy] = VecMoyLogE[indMoy]/Nit;
        VecMoyLogT[indMoy] = VecMoyLogT[indMoy]/Nit;
        VecMoyE[indMoy] = VecMoyE[indMoy]/Nit;
        VecMoyLog0[indMoy] = VecMoyLog0[indMoy]/Nit;
        VecMoyLog1[indMoy] = VecMoyLog1[indMoy]/Nit;
        VecMoyS0[indMoy] = VecMoyS0[indMoy]/Nit;
        VecMoyS1[indMoy] = VecMoyS1[indMoy]/Nit;
        VecMoyXmale[indMoy] = VecMoyXmale[indMoy]/Nit;
        VecMoyY[indMoy] = VecMoyY[indMoy]/Nit;
        VecMoyXfem1[indMoy] = VecMoyXfem1[indMoy]/Nit;
        VecMoyXfem2[indMoy] = VecMoyXfem2[indMoy]/Nit;
        VecMoyH[indMoy] = VecMoyH[indMoy]/Nit;
        //cout << VecMoyXmale[indMoy] << " " << VecMoyY[indMoy] << " " << VecMoyXfem1[indMoy] << " " << VecMoyXfem2[indMoy] << endl;
    }

    vector<vector<double> > DistribH,DistribS;
    int indEchant(0);
    for(indEchant = 0 ; indEchant < Ngen/(100*PasEchant)+1 ; ++indEchant)
    {
        vector<double> DistribHGen(50), DistribSGen(50);
        int indInitVec(0);
        for(indInitVec = 0 ; indInitVec < 50 ; ++indInitVec)
        {
            DistribHGen[indInitVec] = 0;
            DistribSGen[indInitVec] = 0;
        }
        int indDistribIt(0);
        for(indDistribIt = 0 ; indDistribIt < Nit ; ++indDistribIt)
        {
            //cout << log10(VecS1It[indDistribIt][100*indEchant]) << endl;
            double indDistrib(0);
            for(indDistrib = 0 ; indDistrib < 50 ; ++indDistrib)
            {
                if((VecHIt[indDistribIt][100*indEchant] < (indDistrib+1.)/50.) && (VecHIt[indDistribIt][100*indEchant] >= indDistrib/50.))
                {
                    DistribHGen[indDistrib] += 1;
                    //cout << DistribHGen[indDistrib] << endl;
                    //cout << "H+1" << endl;
                }
                if((log10(VecS1It[indDistribIt][100*indEchant]) >= -5.*(indDistrib+1.)/50.) && (log10(VecS1It[indDistribIt][100*indEchant]) < -5.*indDistrib/50.))
                {
                    DistribSGen[indDistrib] += 1;
                    //cout << DistribSGen[indDistrib] << endl;
                    //cout << "S+1" << endl;
                }
                if(VecS1It[indDistribIt][100*indEchant] == 0)
                {
                    DistribSGen[49] += 1;
                }

                //cout << "Distrib " << indDistrib << endl;
                //cout << VecHIt[indDistribIt][100*indEchant] << " " << VecLogS1It[indDistribIt][100*indEchant] << endl;
            }
        }
        DistribH.push_back(DistribHGen);
        DistribS.push_back(DistribSGen);
    }

    /*int indMoyH(0);
    for(indMoyH = 0 ; indMoyH <= (Ngen+TpsDesinit)/PasEchant ; ++indMoyH)
    {
        VecMoyH[indMoyH] = 0;
        VecMoyS[indMoyH] = 0;
        int indTaMereH(0);
        for(indTaMereH = 0 ; indTaMereH < Nit ; ++indTaMereH)
        {
            VecMoyH[indMoyH] = VecMoyH[indMoyH] + VecHIt[indTaMereH][indMoyH];
            cout << VecHIt[indTaMereH][indMoyH] << endl;
            if(SUIVIFIT == "Oui")
            {
                VecMoyW[indMoyH] = VecMoyW[indMoyH] + VecWIt[indTaMereH][indMoyH];
            }
        }
        VecMoyH[indMoyH] = VecMoyH[indMoyH]/Nit;
        VecMoyW[indMoyH] = VecMoyW[indMoyH]/Nit;
    }*/

    int indVecX(0);
    vector<double> X;
    for(indVecX = 0 ; indVecX <= Ngen ; indVecX += PasEchant)
    {
        X.push_back(indVecX);                                                                   //Création du vecteur abscisse de la régression linéaire
    }


    double Pente;
    if(Log == "Oui")
    {
        Pente = LinearRegression(X,VecMoyLogE);
    }
    else if(Log == "Non")
    {
        Pente = LinearRegression(X,VecMoyE);
    }

    if(Results)
    {
        //cout << "Entree des parametres initiaux" << endl;
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
        if(Log == "Oui")
        {
            Results << "Moyenne des logs" << endl;
        }
        else if(Log == "Non")
        {
            Results << "Log de la moyenne" << endl;
        }
        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        Results << "Pente des régressions linéaires du log des forces des Cis-acting factors en fonction du temps : " << Pente << endl;
        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        Results << "Moyennes des traits" << endl;
        Results << endl;
        Results << "T E(X) E(Y) S(X) S(Y) Xmale Y Xfem1 Xfem2 hY" << endl;
        int indH(0);
        for(indH = 0 ; indH < Ngen/PasEchant + 1 ; ++indH)
        {
            Results << indH*PasEchant << " " << VecMoyLogT[indH] << " " << VecMoyLog0[indH] << " " << VecMoyLog1[indH] << " " << VecMoyS0[indH] << " " << VecMoyS1[indH] << " " << VecMoyXmale[indH] << " " << VecMoyY[indH] << " " << VecMoyXfem1[indH] << " " << VecMoyXfem2[indH] << " " << VecMoyH[indH] << endl;
        }
        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        Results << "Evolution des distributions de Hy et Sy" << endl;
        Results << endl;
        int indAffichEchant(0);
        for(indAffichEchant = 0 ; indAffichEchant < Ngen/(100*PasEchant) + 1; ++indAffichEchant)
        {
            Results << "Génération : " << indAffichEchant*100*PasEchant << endl;
            int indAffichDistrib(0);
            for(indAffichDistrib = 0 ; indAffichDistrib < 50 ; ++indAffichDistrib)
            {
                Results << indAffichDistrib*0.02 << " " << DistribH[indAffichEchant][indAffichDistrib] << " " << DistribS[indAffichEchant][indAffichDistrib] << endl;
            }
        }


        if(SUIVIFIT == "Non")
        {
            Results << "Moyennes dans les populations du log des forces des Cis-acting factors : " << endl;

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
                    Results << indGeneration*PasEchant << " " << VecLogEIt[indIteration][indGeneration] << endl;
                }
            }
        }
        else if(SUIVIFIT == "Oui")
        {
            Results << "Moyennes dans les populations du log des forces des Cis-acting factors et des fitness : " << endl;

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
                    Results << indGeneration*PasEchant << " " << VecLogEIt[indIteration][indGeneration] << " " << VecWIt[indIteration][indGeneration] << endl;
                }
            }
        }
    }
}
