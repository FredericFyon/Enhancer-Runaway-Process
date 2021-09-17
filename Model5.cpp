#include "Model5.h"
#include "CycleModel5.h"
#include "MoyFitModel5.h"
#include "Matrice.h"
#include "Statistics.h"
#include "LinearRegression.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>

using namespace std;

void Model5(double Ei, double mutD, double mutE, double mutA, double SigD, double SigE, double h, double S, double Rij, double Rjk, double Self, double I, int Npop, int Ngen1, int Ngen2, int Nit, int TpsDesinit, int PasEchant, string FIT, string ALLELES, string SUIVIFIT, string OutputPathway, string HOMOZYGOTE)
{
    size_t TailleLineEdit = OutputPathway.size() + 1;
    char* FilePathway = new char[TailleLineEdit];
    strncpy(FilePathway,OutputPathway.c_str(),TailleLineEdit);                                  //Construction du fichier de sortie
    ofstream Results(FilePathway,ios::app);

    vector<double> DChrom1(Npop),EChrom1(Npop),SChrom1(Npop),DChrom2(Npop),EChrom2(Npop),SChrom2(Npop);
    vector<double> Empty;

    int indInitPop(0);
    for(indInitPop = 0 ; indInitPop < Npop ; ++indInitPop)
    {
        DChrom1[indInitPop] = 1.0;
        EChrom1[indInitPop] = Ei;
        SChrom1[indInitPop] = 0.;
        DChrom2[indInitPop] = 1.0;
        EChrom2[indInitPop] = Ei;
        SChrom2[indInitPop] = 0.;
    }
    Matrice PopInit(DChrom1,EChrom1,SChrom1,DChrom2,EChrom2,SChrom2,Empty,Empty,Empty,Empty);

    vector<vector<double> > VecLogDIt,VecLogEIt,VecWIt,VecHIt,VecHomIt;

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
            //cout << "Desinitialisation" << endl;
            if(indDesinit == 0)
            {
                Pop = PopInit;
            }
            Pop = CycleModel5(Pop,Ei,0.,0.,mutA,SigD,SigE,h,S,Rij,Rjk,Self,I,Npop,FIT,ALLELES);
            //cout << "Desinit : " << indDesinit << endl;
        }

        vector<double> VecLogD,VecLogE,VecW,VecH,VecHom;

        int indGen(0);
        for(indGen = 0 ; indGen < Ngen1 ; ++indGen)
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
                vector<double> PopH = Pop.getVec7();

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

                //Statistic PopLogE(PopEChrom1,PopEChrom2);
                Statistic PopLogE(PopLogEChrom1,PopLogEChrom2);
                double PopLogEMoy = PopLogE.getMean();
                VecLogE.push_back(PopLogEMoy);

                Statistic PopHMOY(PopH,PopH);
                double PopHMoy = PopHMOY.getMean();
                VecH.push_back(PopHMoy);

                if(SUIVIFIT == "Oui")
                {
                    double WMoy = MoyFitModel5(Ei,PopDChrom1,PopEChrom1,PopSChrom1,PopDChrom2,PopEChrom2,PopSChrom2,h,Npop,FIT,I);
                    VecW.push_back(WMoy);
                }

                double Het(0);
                int indHet(0);
                for(indHet = 0 ; indHet < Npop ; ++indHet)
                {
                    if((PopSChrom1[indHet] != PopSChrom2[indHet]) && (PopEChrom1[indHet] != PopEChrom2[indHet]))
                    {
                        Het = Het + 1;
                    }
                }
                Het = Het/Npop;
                VecHom.push_back(Het);
            }

            Pop = CycleModel5(Pop,Ei,0.,mutE,mutA,SigD,SigE,h,S,Rij,Rjk,Self,I,Npop,FIT,ALLELES);
        }

        for(indGen = Ngen1 ; indGen < Ngen1+Ngen2 ; ++indGen)
        {
            if(indGen%PasEchant == 0)
            {
                vector<double> PopDChrom1 = Pop.getVec1();
                vector<double> PopEChrom1 = Pop.getVec2();
                vector<double> PopSChrom1 = Pop.getVec3();
                vector<double> PopDChrom2 = Pop.getVec4();
                vector<double> PopEChrom2 = Pop.getVec5();
                vector<double> PopSChrom2 = Pop.getVec6();
                vector<double> PopH = Pop.getVec7();

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

                //Statistic PopLogE(PopEChrom1,PopEChrom2);
                Statistic PopLogE(PopLogEChrom1,PopLogEChrom2);
                double PopLogEMoy = PopLogE.getMean();
                VecLogE.push_back(PopLogEMoy);

                Statistic PopHMOY(PopH,PopH);
                double PopHMoy = PopHMOY.getMean();
                VecH.push_back(PopHMoy);

                if(SUIVIFIT == "Oui")
                {
                    double WMoy = MoyFitModel5(Ei,PopDChrom1,PopEChrom1,PopSChrom1,PopDChrom2,PopEChrom2,PopSChrom2,h,Npop,FIT,I);
                    VecW.push_back(WMoy);
                }

                double Het(0);
                int indHet(0);
                for(indHet = 0 ; indHet < Npop ; ++indHet)
                {
                    if((PopSChrom1[indHet] != PopSChrom2[indHet]) && (PopEChrom1[indHet] != PopEChrom2[indHet]))
                    {
                        Het = Het + 1;
                    }
                }
                Het = Het/Npop;
                VecHom.push_back(Het);
            }

            Pop = CycleModel5(Pop,Ei,mutD,mutE,mutA,SigD,SigE,h,S,Rij,Rjk,Self,I,Npop,FIT,ALLELES);
        }

        VecLogDIt.push_back(VecLogD);
        VecLogEIt.push_back(VecLogE);
        VecWIt.push_back(VecW);
        VecHIt.push_back(VecH);
        VecHomIt.push_back(VecHom);
    }

    vector<double> EvolLogE,EvolLogD,EvolW,EvolH,EvolHom;

    if(SUIVIFIT == "Oui")
    {
        int indMoy(0);
        for(indMoy = 0 ; indMoy < VecLogDIt[0].size() ; ++indMoy)
        {
            vector<double> TamponLogE,TamponLogD,TamponW,TamponH,TamponHom;
            int indTampon(0);
            for(indTampon ; indTampon < Nit ; ++indTampon)
            {
                TamponLogE.push_back(VecLogEIt[indTampon][indMoy]);
                TamponLogD.push_back(VecLogDIt[indTampon][indMoy]);
                TamponW.push_back(VecWIt[indTampon][indMoy]);
                TamponH.push_back(VecHIt[indTampon][indMoy]);
                TamponHom.push_back(VecHomIt[indTampon][indMoy]);
            }

            Statistic TAMPONLOGE(TamponLogE,TamponLogE);
            EvolLogE.push_back(TAMPONLOGE.getMean());
            Statistic TAMPONLOGD(TamponLogD,TamponLogD);
            EvolLogD.push_back(TAMPONLOGD.getMean());
            Statistic TAMPONW(TamponW,TamponW);
            EvolW.push_back(TAMPONW.getMean());
            Statistic TAMPONH(TamponH,TamponH);
            EvolH.push_back(TAMPONH.getMean());
            Statistic TAMPONHOM(TamponHom,TamponHom);
            EvolHom.push_back(TAMPONHOM.getMean());
        }
    }

    else if(SUIVIFIT == "Non")
    {
        int indMoy(0);
        for(indMoy = 0 ; indMoy < VecLogDIt[0].size() ; ++indMoy)
        {
            vector<double> TamponLogE,TamponLogD,TamponW,TamponH,TamponHom;
            int indTampon(0);
            for(indTampon ; indTampon < Nit ; ++indTampon)
            {
                TamponLogE.push_back(VecLogEIt[indTampon][indMoy]);
                TamponLogD.push_back(VecLogDIt[indTampon][indMoy]);
                TamponH.push_back(VecHIt[indTampon][indMoy]);
                TamponHom.push_back(VecHomIt[indTampon][indMoy]);
            }

            Statistic TAMPONLOGE(TamponLogE,TamponLogE);
            EvolLogE.push_back(TAMPONLOGE.getMean());
            Statistic TAMPONLOGD(TamponLogD,TamponLogD);
            EvolLogD.push_back(TAMPONLOGD.getMean());
            Statistic TAMPONH(TamponH,TamponH);
            EvolH.push_back(TAMPONH.getMean());
            Statistic TAMPONHOM(TamponHom,TamponHom);
            EvolHom.push_back(TAMPONHOM.getMean());
        }
    }

    if(Results)
    {
        Results << std::endl;
        Results << "_____________________________________________________________________________________________________________" << std::endl;
        Results << std::endl;
        Results << "RESULTATS" << std::endl;
        Results << "_____________________________________________________________________________________________________________" << std::endl;
        Results << std::endl;
        Results << "Valeurs des paramètres : " << std::endl;
        Results << std::endl;
        Results << "Force initiale des Cis-acting factors : " << Ei << std::endl;
        Results << "Coefficient de dominance : " << h << std::endl;
        Results << "Intensité de la sélection au Locus A : " << S << std::endl;
        Results << "Taux de mutation au Locus A : " << mutA << std::endl;
        Results << "Taux de mutation au Locus E : " << mutE << std::endl;
        Results << "Taux de mutation au Locus D : " << mutD << std::endl;
        Results << "Ecart-type de la loi normale régissant les mutations multiplicatives au Locus E : " << SigE << std::endl;
        Results << "Ecart-type de la loi normale régissant les mutations multiplicatives au Locus D : " << SigD << std::endl;
        Results << "Taux de recombinaison entre les locus T et E : " << Rij << std::endl;
        Results << "Taux de recombinaison initial entre les locus E et A : " << Rjk << std::endl;
        Results << "Taux d'autofécondation : " << Self << std::endl;
        Results << "Intensité de la contrainte évolutive : " << I << std::endl;
        Results << "Taille de la population d'individus diploïdes : " << Npop << std::endl;
        Results << "Nombre de générations simulées Phase 1 : " << Ngen1 << std::endl;
        Results << "Nombre de générations simulées Phase 2 : " << Ngen2 << std::endl;
        Results << "Nombre d'itérations du processus évolutif simulées : " << Nit << std::endl;
        Results << "Nombre de générations simulées pour désinitialiser la population : " << TpsDesinit << std::endl;
        Results << "Pas d'échantillonage des valeurs génomiques : " << PasEchant << std::endl;
        if(FIT == "Selection")
        {
            Results << "Sélection active" << std::endl;
        }
        else if(FIT == "Derive")
        {
            Results << "Dérive" << std::endl;
        }
        if(ALLELES == "2")
        {
            Results << "2 allèles au Locus A" << std::endl;
        }
        else if(ALLELES == "Infinite")
        {
            Results << "Nombre infini d'allèles au locus A" << std::endl;
        }
        Results << std::endl;
        Results << "_____________________________________________________________________________________________________________" << std::endl;
        Results << std::endl;
        Results << "Evolution moyenne : " << std::endl;
        Results << std::endl;

        if(SUIVIFIT == "Non")
        {
            int indRes(0);
            if(HOMOZYGOTE == "Oui")
            {
                for(indRes = 0 ; indRes < EvolLogD.size() ; ++indRes)
                {
                    Results << indRes*PasEchant << " " << EvolLogE[indRes] << " " << EvolLogD[indRes] << " " << EvolH[indRes] << " " << EvolHom[indRes] << endl;
                }
            }
            else
            {
                for(indRes = 0 ; indRes < EvolLogD.size() ; ++indRes)
                {
                    Results << indRes*PasEchant << " " << EvolLogE[indRes] << " " << EvolLogD[indRes] << " " << EvolH[indRes] << endl;
                }
            }

            Results << endl;
        }

        else if(SUIVIFIT == "Oui")
        {
            int indRes(0);
            if(HOMOZYGOTE == "Oui")
            {
                for(indRes = 0 ; indRes < EvolLogD.size() ; ++indRes)
                {
                    Results << indRes*PasEchant << " " << EvolLogE[indRes] << " " << EvolLogD[indRes] << " " << EvolH[indRes] << " " << EvolHom[indRes] << " " << EvolW[indRes] << endl;
                }
            }
            else
            {
                for(indRes = 0 ; indRes < EvolLogD.size() ; ++indRes)
                {
                    Results << indRes*PasEchant << " " << EvolLogE[indRes] << " " << EvolLogD[indRes] << " " << EvolH[indRes] << " " << EvolW[indRes] << endl;
                }
            }

            Results << endl;
        }

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
                    Results << indGeneration*PasEchant << " " << VecLogEIt[indIteration][indGeneration]  << " " << VecLogDIt[indIteration][indGeneration] << endl;
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
                    Results << indGeneration*PasEchant << " " << VecLogEIt[indIteration][indGeneration]  << " " << VecLogDIt[indIteration][indGeneration] << " " << VecWIt[indIteration][indGeneration] << endl;
                }
            }
        }
    }
}
