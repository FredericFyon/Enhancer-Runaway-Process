#include "Model3.h"
#include "CycleModel3.h"
#include "MoyFitModel3.h"
#include "Matrice.h"
#include "Statistics.h"
#include "LinearRegression.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include "CalculateVariance.h"
#include "CalculateMean.h"

using namespace std;

void Model3(double Ei, double mutT, double mutE, double mutA, double SigT, double SigE, double h, double S, double Rij, double Rjk, double Self, double I, double Concaveness, int Npop, int Ngen, int Nit, int TpsDesinit, int PasEchant, string FIT, string ALLELES, string SUIVIFIT, string OutputPathway, string SHAPE)
{
    size_t TailleLineEdit = OutputPathway.size() + 1;
    char* FilePathway = new char[TailleLineEdit];
    strncpy(FilePathway,OutputPathway.c_str(),TailleLineEdit);                                  //Construction du fichier de sortie
    ofstream Results(FilePathway,ios::app);

    vector<double> TChrom1(Npop),EChrom1(Npop),SChrom1(Npop),TChrom2(Npop),EChrom2(Npop),SChrom2(Npop);
    vector<double> Empty;

    int indInitPop(0);
    for(indInitPop = 0 ; indInitPop < Npop ; ++indInitPop)
    {
        TChrom1[indInitPop] = pow(10,-log10(Ei));
        EChrom1[indInitPop] = Ei;
        SChrom1[indInitPop] = 0.;
        TChrom2[indInitPop] = pow(10,-log10(Ei));
        EChrom2[indInitPop] = Ei;
        SChrom2[indInitPop] = 0.;
    }
    Matrice PopInit(TChrom1,EChrom1,SChrom1,TChrom2,EChrom2,SChrom2,Empty,Empty,Empty,Empty);

    vector<vector<double> > VecLogTIt,VecLogEIt,VecWIt,VarASEVecIt;
    vector<vector<vector<double> > > VecASEIt;

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
            Pop = CycleModel3(Pop,0.,0.,mutA,SigT,SigE,h,S,Rij,Rjk,Self,I,Concaveness,Npop,FIT,ALLELES,SHAPE);
        }

        vector<double> VecLogT,VecLogE,VecW,VarASEVec;
        vector<vector<double> > VecASETps;

        int indGen(0);
        for(indGen = 0 ; indGen <= Ngen ; ++indGen)
        {
            if(indGen%PasEchant == 0)
            {
                vector<double> PopTChrom1 = Pop.getVec1();
                vector<double> PopEChrom1 = Pop.getVec2();
                vector<double> PopSChrom1 = Pop.getVec3();
                vector<double> PopTChrom2 = Pop.getVec4();
                vector<double> PopEChrom2 = Pop.getVec5();
                vector<double> PopSChrom2 = Pop.getVec6();

                vector<double> PopLogTChrom1,PopLogEChrom1,PopLogTChrom2,PopLogEChrom2;
                int indLog(0);
                for(indLog = 0 ; indLog < Npop ; ++indLog)
                {
                    PopLogTChrom1.push_back(log10(PopTChrom1[indLog]));
                    PopLogEChrom1.push_back(log10(PopEChrom1[indLog]));
                    PopLogTChrom2.push_back(log10(PopTChrom2[indLog]));
                    PopLogEChrom2.push_back(log10(PopEChrom2[indLog]));
                }

                Statistic PopLogT(PopLogTChrom1,PopLogTChrom2);
                double PopLogTMoy = PopLogT.getMean();
                VecLogT.push_back(PopLogTMoy);

                Statistic PopLogE(PopLogEChrom1,PopLogEChrom2);
                double PopLogEMoy = PopLogE.getMean();
                VecLogE.push_back(PopLogEMoy);

                if(SUIVIFIT == "Oui")
                {
                    double WMoy = MoyFitModel3(PopTChrom1,PopEChrom1,PopSChrom1,PopTChrom2,PopEChrom2,PopSChrom2,h,Npop,FIT,I,Concaveness,SHAPE);
                    //Results << WMoy << endl;
                    VecW.push_back(WMoy);
                }
            }

            /*if(indGen%10000 == 0)
            {
                vector<double> PopEChrom1 = Pop.getVec2();
                vector<double> PopEChrom2 = Pop.getVec5();
                vector<double> VecASEIndiv;
                vector<double> VecASEHisto(13,0.);

                int indPop(0);
                for(indPop = 0 ; indPop < Npop ; ++indPop)
                {
                    if(PopEChrom1[indPop] >= PopEChrom2[indPop])
                    {
                        VecASEIndiv.push_back(PopEChrom1[indPop]/(PopEChrom1[indPop]+PopEChrom2[indPop]));
                        //cout << PopEChrom1[indPop]/(PopEChrom1[indPop]+PopEChrom2[indPop]) << endl;
                    }
                    else
                    {
                        VecASEIndiv.push_back(PopEChrom2[indPop]/(PopEChrom1[indPop]+PopEChrom2[indPop]));
                        //cout << PopEChrom2[indPop]/(PopEChrom1[indPop]+PopEChrom2[indPop])<< endl;

                    }
                }

                double MASE = CalculateMean(VecASEIndiv);
                double VarASE = CalculateVariance(VecASEIndiv,MASE);
                VarASEVec.push_back(VarASE);

                int indHisto(0);
                for(indHisto = 0 ; indHisto < Npop ; ++indHisto)
                {
                    if((VecASEIndiv[indHisto] >= 0.5) && (VecASEIndiv[indHisto] < 0.51))
                    {
                        VecASEHisto[0] = VecASEHisto[0] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.51) && (VecASEIndiv[indHisto] < 0.52))
                    {
                        VecASEHisto[1] = VecASEHisto[1] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.52) && (VecASEIndiv[indHisto] < 0.53))
                    {
                        VecASEHisto[2] = VecASEHisto[2] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.53) && (VecASEIndiv[indHisto] < 0.56))
                    {
                        VecASEHisto[3] = VecASEHisto[3] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.56) && (VecASEIndiv[indHisto] < 0.59))
                    {
                        VecASEHisto[4] = VecASEHisto[4] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.59) && (VecASEIndiv[indHisto] < 0.62))
                    {
        }
        else
        {
            EChrom1_Inter = EChrom1_Par1;

            double ProbaRecombinaison2 = TirageUnifReal();
            //Test << ProbaRecombinaison2 << endl;
            if(ProbaRecombinaison2 <= Rjk)
            {
                SChrom1_Inter = SChrom2_Par1;
            }
            else
                        VecASEHisto[5] = VecASEHisto[5] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.62) && (VecASEIndiv[indHisto] < 0.65))
                    {
                        VecASEHisto[6] = VecASEHisto[6] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.65) && (VecASEIndiv[indHisto] < 0.68))
                    {
                        VecASEHisto[7] = VecASEHisto[7] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.68) && (VecASEIndiv[indHisto] < 0.71))
                    {
                        VecASEHisto[8] = VecASEHisto[8] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.71) && (VecASEIndiv[indHisto] < 0.74))
                    {
                        VecASEHisto[9] = VecASEHisto[9] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 74) && (VecASEIndiv[indHisto] < 0.77))
                    {
                        VecASEHisto[10] = VecASEHisto[10] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 77) && (VecASEIndiv[indHisto] < 0.8))
                    {
                        VecASEHisto[11] = VecASEHisto[11] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.8) && (VecASEIndiv[indHisto] < 0.83))
                    {
                        VecASEHisto[12] = VecASEHisto[12] + 1;
                    }
                }

                VecASETps.push_back(VecASEHisto);
            }*/

            Pop = CycleModel3(Pop,mutT,mutE,mutA,SigT,SigE,h,S,Rij,Rjk,Self,I,Concaveness,Npop,FIT,ALLELES,SHAPE);
        }

        VecLogTIt.push_back(VecLogT);
        VecLogEIt.push_back(VecLogE);
        VecWIt.push_back(VecW);
        //VecASEIt.push_back(VecASETps);
        //VarASEVecIt.push_back(VarASEVec);
    }

    vector<double> VecMoyLogT(Ngen/PasEchant+1), VecMoyLogE(Ngen/PasEchant+1);
    int indMoy(0);
    for(indMoy = 0 ; indMoy <= Ngen/PasEchant ; ++indMoy)
    {
        VecMoyLogT[indMoy] = 0;
        VecMoyLogE[indMoy] = 0;
        int indTaMere(0);
        for(indTaMere = 0 ; indTaMere < Nit ; ++indTaMere)
        {
            VecMoyLogT[indMoy] = VecMoyLogT[indMoy] + VecLogTIt[indTaMere][indMoy];
            VecMoyLogE[indMoy] = VecMoyLogE[indMoy] + VecLogEIt[indTaMere][indMoy];
        }
        VecMoyLogT[indMoy] = VecMoyLogT[indMoy]/Nit;
        VecMoyLogE[indMoy] = VecMoyLogE[indMoy]/Nit;
    }

    int indVecX(0);
    vector<double> X;
    for(indVecX = 0 ; indVecX <= Ngen ; indVecX += PasEchant)
    {
        X.push_back(indVecX);                                                                   //Création du vecteur abscisse de la régression linéaire
    }

    //int indPentes(0);
    double PentesT,PentesE;
    PentesT = LinearRegression(X,VecMoyLogT);
    PentesE = LinearRegression(X,VecMoyLogE);


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
        Results << "Taux de mutation au Locus T : " << mutT << endl;
        Results << "Ecart-type de la loi normale régissant les mutations multiplicatives au Locus E : " << SigE << endl;
        Results << "Ecart-type de la loi normale régissant les mutations multiplicatives au Locus T : " << SigT << endl;
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
        Results << "Pentes des régressions linéaires du log des forces des Cis-acting factors en fonction du temps : " << endl;
        Results << "LocusT     LocusE" << endl;
        Results.precision(10);
        Results << fixed;
        Results << PentesT << " " << PentesE << endl;
        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        /*Results << "Variances of Allele-Specific Expression" << endl;
        int indVarASEIt(0);
        for(indVarASEIt = 0 ; indVarASEIt < Nit ; ++indVarASEIt)
        {
            Results << endl;
            Results << "Iteration n°" << indVarASEIt+1 << endl;;
            int indVarASE(0);
            for(indVarASE = 0 ; indVarASE <= Ngen/10000 ; ++indVarASE)
            {
                Results << VarASEVecIt[indVarASEIt][indVarASE] << endl;
            }
        }
        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        Results << "Distribution of Allele-Specific Expression" << endl;
        int indASEIt(0);
        for(indASEIt = 0 ; indASEIt < Nit ; ++indASEIt)
        {
            Results << endl;
            Results << "Iteration n°" << indASEIt+1;
            int indASETps(0);
            for(indASETps = 0 ; indASETps <= Ngen/10000 ; ++indASETps)
            {
                Results << endl;
                Results << "Echant n°" << indASETps << endl;
                int indASEHisto(0);
                for(indASEHisto = 0 ; indASEHisto < 13 ; ++indASEHisto)
                {
                    Results.precision(6);
                    Results << fixed;
                    Results << VecASEIt[indASEIt][indASETps][indASEHisto] << endl;
                }
            }
        }
        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;*/
        if(SUIVIFIT == "Non")
        {
            Results << "Moyennes dans les populations du log des forces des Cis-acting factors : " << endl;
            Results << "Générations    LocusT    LocusE" << endl;

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
                    Results << indGeneration*PasEchant << " " << VecLogTIt[indIteration][indGeneration] << " " << VecLogEIt[indIteration][indGeneration] << endl;
                }
            }
        }
        else if(SUIVIFIT == "Oui")
        {
            Results << "Moyennes dans les populations du log des forces des Cis-acting factors et des fitness : " << endl;
            Results << "Générations    LocusT    LocusE   Fitness" << endl;

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
                    Results << indGeneration*PasEchant << " " << VecLogTIt[indIteration][indGeneration] << " " << VecLogEIt[indIteration][indGeneration] << " " << VecWIt[indIteration][indGeneration] << endl;
                }
            }
        }
    }
}
