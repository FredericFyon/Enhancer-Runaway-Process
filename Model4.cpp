#include "Model4.h"
#include "CycleModel4.h"
#include "MoyFitModel4.h"
#include "Matrice.h"
#include "Statistics.h"
#include "LinearRegression.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>

using namespace std;

void Model4(double Ei, double mutR, double mutE, double mutA, double SigR, double SigE, double h, double S, double Rij, double Rjk, double Self, double I, int Npop, int Ngen1, int Ngen2, int Nit, int TpsDesinit, int PasEchant, string FIT, string ALLELES, string SUIVIFIT, string OutputPathway, string HOMOZYGOTE)
{
    size_t TailleLineEdit = OutputPathway.size() + 1;
    char* FilePathway = new char[TailleLineEdit];
    strncpy(FilePathway,OutputPathway.c_str(),TailleLineEdit);                                  //Construction du fichier de sortie
    ofstream Results(FilePathway,ios::app);

    vector<double> RChrom1(Npop),EChrom1(Npop),SChrom1(Npop),RChrom2(Npop),EChrom2(Npop),SChrom2(Npop);
    vector<double> Empty;

    int indInitPop(0);
    for(indInitPop = 0 ; indInitPop < Npop ; ++indInitPop)
    {
        RChrom1[indInitPop] = Rjk;
        EChrom1[indInitPop] = Ei;
        SChrom1[indInitPop] = 0.;
        RChrom2[indInitPop] = Rjk;
        EChrom2[indInitPop] = Ei;
        SChrom2[indInitPop] = 0.;
    }
    Matrice PopInit(RChrom1,EChrom1,SChrom1,RChrom2,EChrom2,SChrom2,Empty,Empty,Empty,Empty);

    vector<vector<double> > VecLogRIt,VecLogEIt,VecWIt,VecHIt,VecHomIt;

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
            Pop = CycleModel4(Pop,Ei,0.,0.,mutA,SigR,SigE,h,S,Rij,Self,I,Npop,FIT,ALLELES);
        }

        vector<double> VecLogR,VecLogE,VecW,VecH,VecHom;

        int indGen(0);
        for(indGen = 0 ; indGen < Ngen1 ; ++indGen)
        {
            if(indGen%PasEchant == 0)
            {
                vector<double> PopRChrom1 = Pop.getVec1();
                vector<double> PopEChrom1 = Pop.getVec2();
                vector<double> PopSChrom1 = Pop.getVec3();
                vector<double> PopRChrom2 = Pop.getVec4();
                vector<double> PopEChrom2 = Pop.getVec5();
                vector<double> PopSChrom2 = Pop.getVec6();
                vector<double> PopH = Pop.getVec7();

                vector<double> PopLogEChrom1,PopLogEChrom2,PopLogRChrom1,PopLogRChrom2;
                int indLog(0);
                for(indLog = 0 ; indLog < Npop ; ++indLog)
                {
                    PopLogEChrom1.push_back(log10(PopEChrom1[indLog]));
                    PopLogEChrom2.push_back(log10(PopEChrom2[indLog]));
                    PopLogRChrom1.push_back(log10(PopRChrom1[indLog]));
                    PopLogRChrom2.push_back(log10(PopRChrom2[indLog]));
                }

                Statistic PopLogR(PopLogRChrom1,PopLogRChrom2);
                double PopLogRMoy = PopLogR.getMean();
                VecLogR.push_back(PopLogRMoy);

                Statistic PopLogE(PopLogEChrom1,PopLogEChrom2);
                double PopLogEMoy = PopLogE.getMean();
                VecLogE.push_back(PopLogEMoy);

                Statistic POPH(PopH,PopH);
                double PopHMoy = POPH.getMean();
                VecH.push_back(PopHMoy);

                if(SUIVIFIT == "Oui")
                {
                    double WMoy = MoyFitModel4(Ei,PopEChrom1,PopSChrom1,PopEChrom2,PopSChrom2,h,Npop,FIT,I);
                    VecW.push_back(WMoy);
                }

                double Hom(0);
                int indHom(0);
                for(indHom = 0 ; indHom < Npop ; ++indHom)
                {
                    if((PopSChrom1[indHom] != PopSChrom2[indHom]) && (PopEChrom1[indHom] != PopEChrom2[indHom]))
                    {
                        Hom = Hom + 1;
                    }
                }
                Hom = Hom/Npop;
                VecHom.push_back(Hom);
            }

            Pop = CycleModel4(Pop,Ei,0.,mutE,mutA,SigR,SigE,h,S,Rij,Self,I,Npop,FIT,ALLELES);
        }

        for(indGen = Ngen1 ; indGen < Ngen1+Ngen2 ; ++indGen)
        {
            if(indGen%PasEchant == 0)
            {
                vector<double> PopRChrom1 = Pop.getVec1();
                vector<double> PopEChrom1 = Pop.getVec2();
                vector<double> PopSChrom1 = Pop.getVec3();
                vector<double> PopRChrom2 = Pop.getVec4();
                vector<double> PopEChrom2 = Pop.getVec5();
                vector<double> PopSChrom2 = Pop.getVec6();
                vector<double> PopH = Pop.getVec7();

                vector<double> PopLogEChrom1,PopLogEChrom2,PopLogRChrom1,PopLogRChrom2;
                int indLog(0);
                for(indLog = 0 ; indLog < Npop ; ++indLog)
                {
                    PopLogEChrom1.push_back(log10(PopEChrom1[indLog]));
                    PopLogEChrom2.push_back(log10(PopEChrom2[indLog]));
                    PopLogRChrom1.push_back(log10(PopRChrom1[indLog]));
                    PopLogRChrom2.push_back(log10(PopRChrom2[indLog]));
                }

                Statistic PopLogR(PopLogRChrom1,PopLogRChrom2);
                double PopLogRMoy = PopLogR.getMean();
                VecLogR.push_back(PopLogRMoy);

                Statistic PopLogE(PopLogEChrom1,PopLogEChrom2);
                double PopLogEMoy = PopLogE.getMean();
                VecLogE.push_back(PopLogEMoy);

                Statistic POPH(PopH,PopH);
                double PopHMoy = POPH.getMean();
                VecH.push_back(PopHMoy);

                if(SUIVIFIT == "Oui")
                {
                    double WMoy = MoyFitModel4(Ei,PopEChrom1,PopSChrom1,PopEChrom2,PopSChrom2,h,Npop,FIT,I);
                    VecW.push_back(WMoy);
                }

                double Hom(0);
                int indHom(0);
                for(indHom = 0 ; indHom < Npop ; ++indHom)
                {
                    if((PopSChrom1[indHom] != PopSChrom2[indHom]) && (PopEChrom1[indHom] != PopEChrom2[indHom]))
                    {
                        Hom = Hom + 1;
                    }
                }
                Hom = Hom/Npop;
                VecHom.push_back(Hom);
            }

            Pop = CycleModel4(Pop,Ei,mutR,mutE,mutA,SigR,SigE,h,S,Rij,Self,I,Npop,FIT,ALLELES);
        }

        VecLogRIt.push_back(VecLogR);
        VecLogEIt.push_back(VecLogE);
        VecWIt.push_back(VecW);
        VecHIt.push_back(VecH);
        VecHomIt.push_back(VecHom);
    }

    vector<double> EvolLogE,EvolLogR,EvolH,EvolHom;
    int indMoy(0);
    for(indMoy = 0 ; indMoy < VecLogRIt[1].size() ; ++indMoy)
    {
        vector<double> TamponLogE,TamponLogR,TamponH,TamponHom;
        int indTampon(0);
        for(indTampon ; indTampon < Nit ; ++indTampon)
        {
            TamponLogE.push_back(VecLogEIt[indTampon][indMoy]);
            TamponLogR.push_back(VecLogRIt[indTampon][indMoy]);
            TamponH.push_back(VecHIt[indTampon][indMoy]);
            TamponHom.push_back(VecHomIt[indTampon][indMoy]);
        }

        Statistic TAMPONLOGE(TamponLogE,TamponLogE);
        EvolLogE.push_back(TAMPONLOGE.getMean());
        Statistic TAMPONLOGR(TamponLogR,TamponLogR);
        EvolLogR.push_back(TAMPONLOGR.getMean());
        Statistic TAMPONH(TamponH,TamponH);
        EvolH.push_back(TAMPONH.getMean());
        Statistic TAMPONHOM(TamponHom,TamponHom);
        EvolHom.push_back(TAMPONHOM.getMean());
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
        Results << "Taux de mutation au Locus E : " << mutE << endl;
        Results << "Taux de mutation au Locus R : " << mutR << endl;
        Results << "Ecart-type de la loi normale régissant les mutations multiplicatives au Locus E : " << SigE << endl;
        Results << "Ecart-type de la loi normale régissant les mutations multiplicatives au Locus R : " << SigR << endl;
        Results << "Taux de recombinaison entre les locus T et E : " << Rij << endl;
        Results << "Taux de recombinaison initial entre les locus E et A : " << Rjk << endl;
        Results << "Taux d'autofécondation : " << Self << endl;
        Results << "Intensité de la contrainte évolutive : " << I << endl;
        Results << "Taille de la population d'individus diploïdes : " << Npop << endl;
        Results << "Nombre de générations simulées Phase 1: " << Ngen1 << endl;
        Results << "Nombre de générations simulées Phase 2: " << Ngen2 << endl;
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

        if(HOMOZYGOTE == "Oui")
        {
            int indRes(0);
            for(indRes = 0 ; indRes < EvolLogR.size() ; ++indRes)
            {
                Results << indRes*PasEchant << " " << EvolLogE[indRes] << " " << EvolLogR[indRes] << " " << EvolH[indRes] << " " << EvolHom[indRes] << endl;
            }
        }
        else
        {
            int indRes(0);
            for(indRes = 0 ; indRes < EvolLogR.size() ; ++indRes)
            {
                Results << indRes*PasEchant << " " << EvolLogE[indRes] << " " << EvolLogR[indRes] << " " << EvolH[indRes] << endl;
            }
        }

        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;

        if(SUIVIFIT == "Non")
        {
            Results << "Moyennes dans les populations du log des forces des Cis-acting factors : " << endl;
            Results << "Générations    LocusR    LocusE" << endl;

            int indIteration(0);
            for(indIteration = 0 ; indIteration < Nit ; ++ indIteration)
            {
                Results << endl;
                Results << "Itération n°" << indIteration+1 << endl;

                int indGeneration(0);
                for(indGeneration = 0 ; indGeneration <= (Ngen1+Ngen2)/PasEchant ; ++indGeneration)
                {
                    Results.precision(6);
                    Results << fixed;
                    Results << indGeneration*PasEchant << " " << VecLogRIt[indIteration][indGeneration] << " " << VecLogEIt[indIteration][indGeneration] << endl;
                }
            }
        }
        else if(SUIVIFIT == "Oui")
        {
            Results << "Moyennes dans les populations du log des forces des Cis-acting factors et des fitness : " << endl;
            Results << "Générations    LocusR    LocusE   Fitness" << endl;

            int indIteration(0);
            for(indIteration = 0 ; indIteration < Nit ; ++ indIteration)
            {
                Results << endl;
                Results << "Itération n°" << indIteration+1 << endl;

                int indGeneration(0);
                for(indGeneration = 0 ; indGeneration <= (Ngen1+Ngen2)/PasEchant ; ++indGeneration)
                {
                    Results.precision(6);
                    Results << fixed;
                    Results << indGeneration*PasEchant << " " << VecLogRIt[indIteration][indGeneration] << " " << VecLogEIt[indIteration][indGeneration] << " " << VecWIt[indIteration][indGeneration] << endl;
                }
            }
        }
    }
}
