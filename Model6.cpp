#include "Model6.h"
#include "CycleModel6.h"
#include "MoyFitModel6.h"
#include "Matrice.h"
#include "Statistics.h"
#include "LinearRegression.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>

using namespace std;

void Model6(double Ei, double mutD, double mutE, double mutA, double SigD, double SigE, double h, double S, double Rij, double Rjk, double Self, double I, int Npop, int Ngen, int Ngen2, int Nit, int TpsDesinit, int PasEchant, string FIT, string ALLELES, string SUIVIFIT, string OutputPathway, string Modifier, string Log)
{
    //cout << Modifier << endl;

    size_t TailleLineEdit = OutputPathway.size() + 1;
    char* FilePathway = new char[TailleLineEdit];
    strncpy(FilePathway,OutputPathway.c_str(),TailleLineEdit);                                  //Construction du fichier de sortie
    ofstream Results(FilePathway,ios::app);

    vector<double> DChrom1(Npop),EChrom1(Npop),SChrom1(Npop),DChrom2(Npop),EChrom2(Npop),SChrom2(Npop);
    vector<double> Empty;

    if(Modifier == "QProt")
    {
        int indInitPop(0);
        for(indInitPop = 0 ; indInitPop < Npop ; ++indInitPop)
        {
            DChrom1[indInitPop] = 1.00;
            EChrom1[indInitPop] = Ei;
            SChrom1[indInitPop] = 0.;
            DChrom2[indInitPop] = 1.00;
            EChrom2[indInitPop] = Ei;
            SChrom2[indInitPop] = 0.;
        }
    }
    else if(Modifier == "ExcesQProt")
    {
        int indInitPop(0);
        for(indInitPop = 0 ; indInitPop < Npop ; ++indInitPop)
        {
            DChrom1[indInitPop] = 1.00;
            EChrom1[indInitPop] = Ei;
            SChrom1[indInitPop] = 0.;
            DChrom2[indInitPop] = 1.00;
            EChrom2[indInitPop] = Ei;
            SChrom2[indInitPop] = 0.;
        }
    }

    Matrice PopInit(DChrom1,EChrom1,SChrom1,DChrom2,EChrom2,SChrom2,Empty,Empty,Empty,Empty);

    vector<vector<double> > VecLogDIt,VecLogEIt,VecWIt,VecDIt,VecEIt;

    double indIt(0);
    for(indIt = 0 ; indIt < Nit ; ++indIt)
    {
        if(Results)
        {
            Results << "Iteration n°" << indIt+1 << " in progress..." << endl;
        }
        //cout << indIt+1 << endl;

        Matrice Pop(Empty,Empty,Empty,Empty,Empty,Empty,Empty,Empty,Empty,Empty);

        int indDesinit(0);
        for(indDesinit = 0 ; indDesinit < TpsDesinit ; ++indDesinit)
        {
            if(indDesinit == 0)
            {
                Pop = PopInit;
            }
            Pop = CycleModel6(Pop,Ei,0.,0.,mutA,SigD,SigE,h,S,Rij,Rjk,Self,I,Npop,FIT,ALLELES,Modifier);
        }

        vector<double> VecLogD,VecLogE,VecW,VecD,VecE;

        int indGen(0);
        for(indGen = 0 ; indGen < Ngen ; ++indGen)
        {
            if(indGen%PasEchant == 0)
            {
                //cout << indGen << endl;
                vector<double> PopDChrom1 = Pop.getVec1();
                vector<double> PopEChrom1 = Pop.getVec2();
                vector<double> PopSChrom1 = Pop.getVec3();
                vector<double> PopDChrom2 = Pop.getVec4();
                vector<double> PopEChrom2 = Pop.getVec5();
                vector<double> PopSChrom2 = Pop.getVec6();

                vector<double> PopLogEChrom1,PopLogEChrom2,PopLogDChrom1,PopLogDChrom2;
                int indLog(0);
                for(indLog = 0 ; indLog < Npop ; ++indLog)
                {
                    PopLogEChrom1.push_back(log10(PopEChrom1[indLog]));
                    PopLogEChrom2.push_back(log10(PopEChrom2[indLog]));
                    PopLogDChrom1.push_back(log10(PopDChrom1[indLog]));
                    PopLogDChrom2.push_back(log10(PopDChrom2[indLog]));
                }

                //Statistic PopD(PopDChrom1,PopDChrom2);
                Statistic PopLogD(PopLogDChrom1,PopLogDChrom2);
                double PopLogDMoy = PopLogD.getMean();
                VecLogD.push_back(PopLogDMoy);
                Statistic PopD(PopDChrom1,PopDChrom2);
                double PopDMoy = PopD.getMean();
                VecD.push_back(PopDMoy);

                Statistic PopE(PopEChrom1,PopEChrom2);
                double PopEMoy = PopE.getMean();
                VecE.push_back(PopEMoy);
                Statistic PopLogE(PopLogEChrom1,PopLogEChrom2);
                double PopLogEMoy = PopLogE.getMean();
                VecLogE.push_back(PopLogEMoy);

                if(SUIVIFIT == "Oui")
                {
                    double WMoy = MoyFitModel6(Ei,PopDChrom1,PopEChrom1,PopSChrom1,PopDChrom2,PopEChrom2,PopSChrom2,h,Npop,FIT,I,Modifier);
                    VecW.push_back(WMoy);
                }
            }

            Pop = CycleModel6(Pop,Ei,0.,mutE,mutA,SigD,SigE,h,S,Rij,Rjk,Self,I,Npop,FIT,ALLELES,Modifier);
        }

        //cout << "Fin initialisation" << endl;

        for(indGen = Ngen ; indGen <= Ngen2+Ngen ; ++indGen)
        {
            if(indGen%PasEchant == 0)
            {
                //cout << indGen << endl;
                vector<double> PopDChrom1 = Pop.getVec1();
                vector<double> PopEChrom1 = Pop.getVec2();
                vector<double> PopSChrom1 = Pop.getVec3();
                vector<double> PopDChrom2 = Pop.getVec4();
                vector<double> PopEChrom2 = Pop.getVec5();
                vector<double> PopSChrom2 = Pop.getVec6();

                vector<double> PopLogEChrom1,PopLogEChrom2,PopLogDChrom1,PopLogDChrom2;
                int indLog(0);
                for(indLog = 0 ; indLog < Npop ; ++indLog)
                {
                    PopLogEChrom1.push_back(log10(PopEChrom1[indLog]));
                    PopLogEChrom2.push_back(log10(PopEChrom2[indLog]));
                    PopLogDChrom1.push_back(log10(PopDChrom1[indLog]));
                    PopLogDChrom2.push_back(log10(PopDChrom2[indLog]));
                }

                Statistic PopLogD(PopLogDChrom1,PopLogDChrom2);
                double PopLogDMoy = PopLogD.getMean();
                VecLogD.push_back(PopLogDMoy);
                Statistic PopD(PopDChrom1,PopDChrom2);
                double PopDMoy = PopD.getMean();
                VecD.push_back(PopDMoy);

                Statistic PopE(PopEChrom1,PopEChrom2);
                double PopEMoy = PopE.getMean();
                VecE.push_back(PopEMoy);
                Statistic PopLogE(PopLogEChrom1,PopLogEChrom2);
                double PopLogEMoy = PopLogE.getMean();
                VecLogE.push_back(PopLogEMoy);

                if(SUIVIFIT == "Oui")
                {
                    double WMoy = MoyFitModel6(Ei,PopDChrom1,PopEChrom1,PopSChrom1,PopDChrom2,PopEChrom2,PopSChrom2,h,Npop,FIT,I,Modifier);
                    VecW.push_back(WMoy);
                }
            }

            Pop = CycleModel6(Pop,Ei,mutD,mutE,mutA,SigD,SigE,h,S,Rij,Rjk,Self,I,Npop,FIT,ALLELES,Modifier);
        }

        VecLogDIt.push_back(VecLogD);
        VecLogEIt.push_back(VecLogE);
        VecDIt.push_back(VecD);
        VecEIt.push_back(VecE);
        VecWIt.push_back(VecW);
    }

    //cout << "Fin cycles de vie" << endl;

    vector<double> EvolLogE,EvolLogD,EvolW,EvolE,EvolD;

    if(SUIVIFIT == "Oui")
    {
        int indMoy(0);
        for(indMoy = 0 ; indMoy < VecLogDIt[1].size() ; ++indMoy)
        {
            vector<double> TamponLogE,TamponLogD,TamponW,TamponD,TamponE;
            int indTampon(0);
            for(indTampon ; indTampon < Nit ; ++indTampon)
            {
                TamponLogE.push_back(VecLogEIt[indTampon][indMoy]);
                TamponE.push_back(VecEIt[indTampon][indMoy]);
                TamponLogD.push_back(VecLogDIt[indTampon][indMoy]);
                TamponD.push_back(VecDIt[indTampon][indMoy]);
                TamponW.push_back(VecWIt[indTampon][indMoy]);
            }

            Statistic TAMPONLOGE(TamponLogE,TamponLogE);
            EvolLogE.push_back(TAMPONLOGE.getMean());
            Statistic TAMPONE(TamponE,TamponE);
            EvolE.push_back(log10(TAMPONE.getMean()));
            Statistic TAMPONLOGD(TamponLogD,TamponLogD);
            EvolLogD.push_back(TAMPONLOGD.getMean());
            Statistic TAMPOND(TamponD,TamponD);
            EvolD.push_back(log10(TAMPOND.getMean()));
            Statistic TAMPONW(TamponW,TamponW);
            EvolW.push_back(TAMPONW.getMean());
        }
    }

    else if(SUIVIFIT == "Non")
    {
        int indMoy(0);
        for(indMoy = 0 ; indMoy < (Ngen2+Ngen)/PasEchant + 1 ; ++indMoy)
        {
            //cout << "Réorganisation : " << indMoy << endl;
            vector<double> TamponLogE,TamponLogD,TamponW,TamponE,TamponD;
            int indTampon(0);
            for(indTampon ; indTampon < Nit ; ++indTampon)
            {
                TamponLogE.push_back(VecLogEIt[indTampon][indMoy]);
                TamponE.push_back(VecEIt[indTampon][indMoy]);
                TamponLogD.push_back(VecLogDIt[indTampon][indMoy]);
                TamponD.push_back(VecDIt[indTampon][indMoy]);
            }

            Statistic TAMPONLOGE(TamponLogE,TamponLogE);
            EvolLogE.push_back(TAMPONLOGE.getMean());
            Statistic TAMPONE(TamponE,TamponE);
            EvolE.push_back(log10(TAMPONE.getMean()));
            Statistic TAMPONLOGD(TamponLogD,TamponLogD);
            EvolLogD.push_back(TAMPONLOGD.getMean());
            Statistic TAMPOND(TamponD,TamponD);
            EvolD.push_back(log10(TAMPOND.getMean()));
        }
    }

    //cout << "Debut ecriture" << endl;

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
        Results << "Taux de mutation au Locus M : " << mutD << endl;
        Results << "Ecart-type de la loi normale régissant les mutations multiplicatives au Locus E : " << SigE << endl;
        Results << "Ecart-type de la loi normale régissant les mutations multiplicatives au Locus D : " << SigD << endl;
        Results << "Taux de recombinaison entre les locus T et E : " << Rij << endl;
        Results << "Taux de recombinaison initial entre les locus E et A : " << Rjk << endl;
        Results << "Taux d'autofécondation : " << Self << endl;
        Results << "Intensité de la contrainte évolutive : " << I << endl;
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
        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        Results << "Evolution moyenne : " << endl;
        Results << endl;

        if(SUIVIFIT == "Non")
        {
            int indRes(0);
            for(indRes = 0 ; indRes < EvolD.size() ; ++indRes)
            {
                if(Log == "Oui")
                {
                    Results << indRes*PasEchant << " " << EvolLogE[indRes] << " " << EvolLogD[indRes] << endl;
                }
                else if(Log == "Non")
                {
                    Results << indRes*PasEchant << " " << EvolE[indRes] << " " << EvolD[indRes] << endl;
                }
            }
            Results << endl;
        }

        else if(SUIVIFIT == "Oui")
        {
            int indRes(0);
            for(indRes = 0 ; indRes < EvolD.size() ; ++indRes)
            {
                if(Log == "Oui")
                {
                    Results << indRes*PasEchant << " " << EvolLogE[indRes] << " " << EvolLogD[indRes] << " " << EvolW[indRes] << endl;
                }
                else if(Log == "Non")
                {
                    Results << indRes*PasEchant << " " << EvolE[indRes] << " " << EvolD[indRes] << " " << EvolW[indRes] << endl;
                }
            }
            Results << endl;
        }

        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;

        if(SUIVIFIT == "Non")
        {
            Results << "Moyennes dans les populations du log des forces des Cis-acting factors : " << endl;
            Results << "Générations    LocusM    LocusE" << endl;

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
                    Results << indGeneration*PasEchant << " " << VecLogEIt[indIteration][indGeneration]  << " " << VecLogDIt[indIteration][indGeneration] << endl;
                }
            }
        }
        else if(SUIVIFIT == "Oui")
        {
            Results << "Moyennes dans les populations du log des forces des Cis-acting factors et des fitness : " << endl;
            Results << "Générations    LocusM    LocusE   Fitness" << endl;

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
                    Results << indGeneration*PasEchant << " " << VecLogEIt[indIteration][indGeneration]  << " " << VecLogDIt[indIteration][indGeneration] << " " << VecWIt[indIteration][indGeneration] << endl;
                }
            }
        }
    }
}
