#include "Model1.h"
#include "CycleModel1.h"
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

using namespace std;

void Model1(double Ei, double mutE, double mutA, double SigE, double h, double S, double Rij, double Rjk, double Self, double Clone, double Automixis, int Npop, int Ngen, int Nit, int TpsDesinit, int PasEchant, string FIT, string ALLELES, string SUIVIFIT, string OutputPathway, string Log, string CANCER, string HOMOZYGOTE, string AUTOMIXIS)
{
    size_t TailleLineEdit = OutputPathway.size() + 1;
    char* FilePathway = new char[TailleLineEdit];
    strncpy(FilePathway,OutputPathway.c_str(),TailleLineEdit);                                  //Construction du fichier de sortie
    ofstream Results(FilePathway,ios::app);

    ofstream Test("/home/frederic/Bureau/GEFSimulator_CB/bin/Debug/Test.txt",ios::app);

    vector<double> EChrom1(Npop),SChrom1(Npop),EChrom2(Npop),SChrom2(Npop);                     //Création des 4 vecteurs représentant les 4 valeurs génomiques de chaque individu
    vector<double> Empty;                                                                       //Création d'un vecteur vide

    int indInitPop(0);
    for(indInitPop = 0 ; indInitPop < Npop ; ++indInitPop)
    {
        EChrom1[indInitPop] = Ei;
        EChrom2[indInitPop] = Ei;
        SChrom1[indInitPop] = 0.;                                                               //Initiation des valeurs génomiques : population monomorphique
        SChrom2[indInitPop] = 0.;
    }
    Matrice PopInit(EChrom1,SChrom1,EChrom2,SChrom2,Empty,Empty,Empty,Empty,Empty,Empty);       //Matrice représentant la population initiale

    vector<vector<double> > VecLogEIt,VecWIt,VarASEVecIt,VecEIt,VecHIt,VecSIt,VecLogEplusIt,VecLogEmoinsIt,VecSplusIt,VecSmoinsIt,VecHomIt;
    vector<vector<vector<double> > > VecASEIt;                                                   //Création des matrices stockant les forces des Cis-acting factors et les fitness au cours
                                                                                                //des générations et pour les Nit itérations
    double indIt(0);
    for(indIt = 0 ; indIt < Nit ; ++indIt)                                                      //Boucle des Nit itérations de l'évolution d'une même population initiale
    {
        vector<double> VecLogE,VecW,VarASEVec,VecE,VecH,VecS,VecLogEplus,VecLogEmoins,VecSplus,VecSmoins,VecHom;
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
            Pop = CycleModel1(Pop,0.,mutA,SigE,h,S,Rij,Rjk,Self,Clone,Automixis,Npop,FIT,ALLELES,AUTOMIXIS);
            /*if(indDesinit%PasEchant == 0)
            {
                vector<double> PopEChrom1 = Pop.getVec1();
                vector<double> PopEChrom2 = Pop.getVec3();
                vector<double> PopLogEChrom1,PopLogEChrom2;
                int indLog(0);
                for(indLog = 0 ; indLog < Npop ; ++indLog)
                {
                    PopLogEChrom1.push_back(log10(PopEChrom1[indLog]));                         //Passage des données au logarithme décimal
                    PopLogEChrom2.push_back(log10(PopEChrom2[indLog]));
                }
                Statistic PopLogE(PopLogEChrom1,PopLogEChrom2);
                double PopLogEMoy = PopLogE.getMean();
                VecLogE.push_back(PopLogEMoy);
                vector<double> PopH = Pop.getVec5();
                Statistic PopHMOY(PopH,PopH);
                double PopHMoy = PopHMOY.getMean();
                VecH.push_back(PopHMoy);
                vector<double> PopS1 = Pop.getVec2();
                vector<double> PopS2 = Pop.getVec4();
                Statistic PopSMOY(PopS1,PopS2);
                double PopSMoy = PopSMOY.getMean();
                VecS.push_back(PopSMoy);
                if(SUIVIFIT == "Oui")
                {
                    double WMoy = MoyFitModel1(Pop.getVec1(),Pop.getVec2(),Pop.getVec3(),Pop.getVec4(),h,Npop,FIT);    //Stockage à intervalles réguliers de la moyenne de la fitness dans la population
                    //Results << WMoy << endl;
                    VecW.push_back(WMoy);
                }
            }*/
        }
                                                                //Création des vecteurs stockant les forces des Cis-acting factors et les fitness au cours
        Test << "Fin Desinit" << endl;                                                                                        //des générations
        int indGen(0);
        for(indGen = 0 ; indGen <= Ngen ; ++indGen)
        {
            /*if(Results)
            {
                int indPop(0);
                vector<double> Repro = Pop.getVec6();
                for(indPop = 0 ; indPop < Npop ; ++indPop)
                {
                    Results << Repro[indPop] << endl;
                }
            }*/
            if(indGen%PasEchant == 0)
            {
                //cout << indGen << endl;
                vector<double> PopEChrom1 = Pop.getVec1();
                vector<double> PopSChrom1 = Pop.getVec2();
                vector<double> PopEChrom2 = Pop.getVec3();                                      //Extraction à intervalle régulier des valeurs génomiques de la population
                vector<double> PopSChrom2 = Pop.getVec4();
                vector<double> PopH = Pop.getVec5();

                vector<double> PopLogEChrom1,PopLogEChrom2;
                int indLog(0);
                for(indLog = 0 ; indLog < Npop ; ++indLog)
                {
                    PopLogEChrom1.push_back(log10(PopEChrom1[indLog]));
                    PopLogEChrom2.push_back(log10(PopEChrom2[indLog]));
                }

                Statistic PopLogE(PopLogEChrom1,PopLogEChrom2);
                double PopLogEMoy = PopLogE.getMean();
                VecLogE.push_back(PopLogEMoy);
                Statistic PopE(PopEChrom1,PopEChrom2);
                double PopEMoy = log10(PopE.getMean());
                VecE.push_back(PopEMoy);
                Statistic PopHMOY(PopH,PopH);
                double PopHMoy = PopHMOY.getMean();
                VecH.push_back(PopHMoy);
                Statistic PopSMOY(PopSChrom1,PopSChrom2);
                double PopSMoy = PopSMOY.getMean();
                VecS.push_back(PopSMoy);

                if(CANCER == "Oui")
                {
                    vector<double> LogEplus;
                    vector<double> LogEmoins;
                    vector<double> Splus;
                    vector<double> Smoins;
                    int indCancer(0);
                    for(indCancer = 0 ; indCancer < Npop ; ++indCancer)
                    {
                        if(PopEChrom1[indCancer] >= PopEChrom2[indCancer])
                        {
                            LogEplus.push_back(log10(PopEChrom1[indCancer]));
                            LogEmoins.push_back(log10(PopEChrom2[indCancer]));
                            Splus.push_back(PopSChrom1[indCancer]);
                            Smoins.push_back(PopSChrom2[indCancer]);
                        }
                        else
                        {
                            LogEplus.push_back(log10(PopEChrom2[indCancer]));
                            LogEmoins.push_back(log10(PopEChrom1[indCancer]));
                            Splus.push_back(PopSChrom2[indCancer]);
                            Smoins.push_back(PopSChrom1[indCancer]);
                        }
                    }

                    Statistic PopLogEPLUS(LogEplus,LogEplus);
                    double PopLogEplus = PopLogEPLUS.getMean();
                    VecLogEplus.push_back(PopLogEplus);

                    Statistic PopLogEMOINS(LogEmoins,LogEmoins);
                    double PopLogEmoins = PopLogEMOINS.getMean();
                    VecLogEmoins.push_back(PopLogEmoins);

                    Statistic PopSPLUS(Splus,Splus);
                    double PopSPlus = PopSPLUS.getMean();
                    VecSplus.push_back(PopSPlus);

                    Statistic PopSMOINS(Smoins,Smoins);
                    double PopSMoins = PopSMOINS.getMean();
                    VecSmoins.push_back(PopSMoins);
                }

                if(HOMOZYGOTE == "Oui")
                {
                    double Hom(0);
                    int indHom(0);
                    for(indHom = 0 ; indHom < Npop ; ++indHom)
                    {
                        if(PopSChrom1[indHom] == PopSChrom2[indHom])
                        {
                            Hom = Hom + 1;
                        }

                        if(PopEChrom1[indHom] == PopEChrom2[indHom])
                        {
                            Hom = Hom + 1;
                        }
                    }
                    Hom = Hom/(2*Npop);
                    //Test  << Hom << endl;
                    VecHom.push_back(Hom);
                }

                if(SUIVIFIT == "Oui")
                {
                    double WMoy = MoyFitModel1(PopEChrom1,PopSChrom1,PopEChrom2,PopSChrom2,h,Npop,FIT);    //Stockage à intervalles réguliers de la moyenne de la fitness dans la population
                    //Test << WMoy << endl;
                    VecW.push_back(WMoy);
                }
            }



          /*  if(indGen%10000 == 0)
            {
                vector<double> PopEChrom1 = Pop.getVec1();
                vector<double> PopEChrom2 = Pop.getVec3();
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
                    if((VecASEIndiv[indHisto] >= 0.5) && (VecASEIndiv[indHisto] < 0.52))
                    {
                        VecASEHisto[0] = VecASEHisto[0] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.52) && (VecASEIndiv[indHisto] < 0.54))
                    {
                        VecASEHisto[1] = VecASEHisto[1] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.54) && (VecASEIndiv[indHisto] < 0.56))
                    {
                        VecASEHisto[2] = VecASEHisto[2] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.56) && (VecASEIndiv[indHisto] < 0.58))
                    {
                        VecASEHisto[3] = VecASEHisto[3] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.58) && (VecASEIndiv[indHisto] < 0.6))
                    {
                        VecASEHisto[4] = VecASEHisto[4] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.6) && (VecASEIndiv[indHisto] < 0.65))
                    {
                        VecASEHisto[5] = VecASEHisto[5] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.65) && (VecASEIndiv[indHisto] < 0.7))
                    {
                        VecASEHisto[6] = VecASEHisto[6] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.7) && (VecASEIndiv[indHisto] < 0.75))
                    {
                        VecASEHisto[7] = VecASEHisto[7] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.75) && (VecASEIndiv[indHisto] < 0.8))
                    {
                        VecASEHisto[8] = VecASEHisto[8] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.8) && (VecASEIndiv[indHisto] < 0.85))
                    {
                        VecASEHisto[9] = VecASEHisto[9] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.85) && (VecASEIndiv[indHisto] < 0.9))
                    {
                        VecASEHisto[10] = VecASEHisto[10] + 1;
                    }
                    else if((VecASEIndiv[indHisto] >= 0.9) && (VecASEIndiv[indHisto] < 0.95))
                    {
                        VecASEHisto[11] = VecASEHisto[11] + 1;
                    }
                    else
                    {
                        VecASEHisto[12] = VecASEHisto[12] + 1;
                    }
                }

                VecASETps.push_back(VecASEHisto);
            }*/

            Pop = CycleModel1(Pop,mutE,mutA,SigE,h,S,Rij,Rjk,Self,Clone,Automixis,Npop,FIT,ALLELES,AUTOMIXIS);
        }

        VecLogEIt.push_back(VecLogE);
        VecLogEplusIt.push_back(VecLogEplus);
        VecLogEmoinsIt.push_back(VecLogEmoins);
        VecSplusIt.push_back(VecSplus);
        VecSmoinsIt.push_back(VecSmoins);
        VecEIt.push_back(VecE);
        VecWIt.push_back(VecW);
        VecHIt.push_back(VecH);
        VecSIt.push_back(VecS);                                                             //Stockage à chaque itération de l'évolution de la moyenne dans la pop du log des forces des
        VecASEIt.push_back(VecASETps);
        VarASEVecIt.push_back(VarASEVec);
        VecHomIt.push_back(VecHom);                                                                                      //Cis-acting factors et de la moyenne dans la pop des fitness
    }

    vector<double> VecMoyLogE((Ngen+TpsDesinit)/PasEchant+1), VecMoyLogEmoins(Ngen/PasEchant+1), VecMoySplus(Ngen/PasEchant+1), VecMoySmoins(Ngen/PasEchant+1), VecMoyLogEplus(Ngen/PasEchant+1), VecMoyE(Ngen/PasEchant+1),VecMoyH((Ngen+TpsDesinit)/PasEchant+1),VecMoyW((Ngen+TpsDesinit)/PasEchant + 1),VecMoyS((Ngen+TpsDesinit)/PasEchant+1),VecMoyHom(Ngen/PasEchant+1);
    int indMoy(0);
    for(indMoy = 0 ; indMoy <= Ngen/PasEchant ; ++indMoy)
    {
        VecMoyE[indMoy] = 0;
        VecMoyLogEplus[indMoy] = 0;
        VecMoyLogEmoins[indMoy] = 0;
        VecMoySplus[indMoy] = 0;
        VecMoySmoins[indMoy] = 0;
        VecMoyHom[indMoy] = 0;
        int indTaMere(0);
        for(indTaMere = 0 ; indTaMere < Nit ; ++indTaMere)
        {
            VecMoyE[indMoy] = VecMoyE[indMoy] + VecEIt[indTaMere][indMoy];
            if(CANCER == "Oui")
            {
                VecMoyLogEplus[indMoy] = VecMoyLogEplus[indMoy] + VecLogEplusIt[indTaMere][indMoy];
                VecMoyLogEmoins[indMoy] = VecMoyLogEmoins[indMoy] + VecLogEmoinsIt[indTaMere][indMoy];
                VecMoySplus[indMoy] = VecMoySplus[indMoy] + VecSplusIt[indTaMere][indMoy];
                VecMoySmoins[indMoy] = VecMoySmoins[indMoy] + VecSmoinsIt[indTaMere][indMoy];
            }

            if(HOMOZYGOTE == "Oui")
            {
                VecMoyHom[indMoy] = VecMoyHom[indMoy] + VecHomIt[indTaMere][indMoy];
            }
        }
        VecMoyE[indMoy] = VecMoyE[indMoy]/Nit;
        VecMoyLogEplus[indMoy] = VecMoyLogEplus[indMoy]/Nit;
        VecMoyLogEmoins[indMoy] = VecMoyLogEmoins[indMoy]/Nit;
        VecMoySplus[indMoy] = VecMoySplus[indMoy]/Nit;
        VecMoySmoins[indMoy] = VecMoySmoins[indMoy]/Nit;
        VecMoyHom[indMoy] = VecMoyHom[indMoy]/Nit;
    }

    /*if(Results)
    {
        Results << "Changements vecteurs" << endl;
    }*/

    int indMoyH(0);
    for(indMoyH = 0 ; indMoyH <= (Ngen+TpsDesinit)/PasEchant ; ++indMoyH)
    {
        VecMoyLogE[indMoyH] = 0;
        VecMoyH[indMoyH] = 0;
        VecMoyW[indMoyH] = 0;
        VecMoyS[indMoyH] = 0;
        int indTaMereH(0);
        for(indTaMereH = 0 ; indTaMereH < Nit ; ++indTaMereH)
        {
            VecMoyH[indMoyH] = VecMoyH[indMoyH] + VecHIt[indTaMereH][indMoyH];
            if(SUIVIFIT == "Oui")
            {
                VecMoyW[indMoyH] = VecMoyW[indMoyH] + VecWIt[indTaMereH][indMoyH];
            }
            VecMoyS[indMoyH] = VecMoyS[indMoyH] + VecSIt[indTaMereH][indMoyH];
            VecMoyLogE[indMoyH] = VecMoyLogE[indMoyH] + VecLogEIt[indTaMereH][indMoyH];
        }
        VecMoyLogE[indMoyH] = VecMoyLogE[indMoyH]/Nit;
        VecMoyH[indMoyH] = VecMoyH[indMoyH]/Nit;
        VecMoyW[indMoyH] = VecMoyW[indMoyH]/Nit;
        VecMoyS[indMoyH] = VecMoyS[indMoyH]/Nit;
    }

    //Results << "Moyennes calculees" << endl;

    int indVecX(0);
    vector<double> X;
    for(indVecX = 0 ; indVecX <= Ngen ; indVecX += PasEchant)
    {
        X.push_back(indVecX);                                                                   //Création du vecteur abscisse de la régression linéaire
    }

    int indPentes(0);
    vector<double> Pentes;
    for(indPentes = 0 ; indPentes < Nit ; ++indPentes)
    {
        //cout << "calcul pente n°" << indPentes+1 << endl;
        if(Log == "Oui")
        {
            Pentes.push_back(LinearRegression(X,VecLogEIt[indPentes]));
        }
        else if(Log == "Non")
        {
            Pentes.push_back(LinearRegression(X,VecEIt[indPentes]));
        }                                                                                    //Calcul des pentes des régressions linéraires du log des forces des Cis-acting factors en fonction du temps
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
        Results << "Fréquence Automixis : " << Automixis << endl;
        Results << "Fréquence Cloning : " << Clone << endl;
        Results << "Fréquence Selfing : " << Self << endl;
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
        /*int indPente(0);
        for(indPente = 0 ; indPente < Nit ; ++indPente)
        {
            //cout << "Entree pentes n°" << indPente+1 << endl;
            Results.precision(10);
            Results << fixed;
            Results << Pentes[indPente] << endl;
        }*/

        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        if(CANCER == "Oui")
        {
            Results << "Asymetrie clonale" << endl;
            Results << "Generation E+ E- S+ S-" << endl;
            int indClone(0);
            for(indClone = 0 ; indClone < (Ngen+TpsDesinit)/PasEchant + 1 ; ++indClone)
            {
                Results << indClone*PasEchant << " " << VecMoyLogEplus[indClone] << " " << VecMoyLogEmoins[indClone] << " " << VecMoySplus[indClone] << " " << VecMoySmoins[indClone] << endl;
            }
            Results << endl;
            Results << "_____________________________________________________________________________________________________________" << endl;
            Results << endl;
        }
        if(HOMOZYGOTE == "Oui")
        {
            Results << "Homozygosite" << endl;
            int indHom(0);
            for(indHom = 0 ; indHom < Ngen/PasEchant + 1 ; ++indHom)
            {
                Results << indHom*PasEchant << " " << VecMoyHom[indHom] << endl;
            }
            Results << endl;
            Results << "_____________________________________________________________________________________________________________" << endl;
            Results << endl;
        }
        Results << "Dominance et fitness moyenne" << endl;
        int indH(0);
        for(indH = 0 ; indH < (Ngen+TpsDesinit)/PasEchant + 1 ; ++indH)
        {
            Results << indH*PasEchant << " " << VecMoyLogE[indH] << " " << VecMoyS[indH] << " " << VecMoyH[indH] << " " << VecMoyW[indH] << endl;
        }
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
        }*/
        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;

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

