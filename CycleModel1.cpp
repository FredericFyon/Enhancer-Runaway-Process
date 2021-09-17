#include "CycleModel1.h"
#include "Matrice.h"
#include "Vecteur.h"
#include "RecombinaisonModel1.h"
#include "RecombinaisonAutomixis.h"
#include "SelectionModel1.h"
#include <string>
#include <fstream>

#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>

#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions.hpp>

using namespace std;

Matrice CycleModel1(Matrice Pop, double mutE, double mutA, double SigE, double h, double S, double Rij, double Rjk, double Self, double Clone, double Automixis, int Npop, string FIT, string ALLELES, string AUTOMIXIS)
{
    //ofstream Test("/home/frederic/Bureau/GEFSimulator_CB/bin/Debug/Test.txt",ios::app);

    vector<double> EChrom1 = Pop.getVec1();
    vector<double> SChrom1 = Pop.getVec2();
    vector<double> EChrom2 = Pop.getVec3();
    vector<double> SChrom2 = Pop.getVec4();
    vector<double> ECHROM1(Npop),ECHROM2(Npop),SCHROM1(Npop),SCHROM2(Npop);
    vector<double> W,H,Empty;

    typedef boost::normal_distribution<float> NormalDistribution;
    typedef boost::exponential_distribution<double> ExponentialDistribution;
    typedef boost::uniform_int<int> UniformDistributionInt;
    typedef boost::uniform_real<double> UniformDistributionReal;
    typedef boost::poisson_distribution<int> PoissonDistribution;

    typedef boost::mt19937 RandomGenerator;

    typedef boost::variate_generator<RandomGenerator&, NormalDistribution> GaussianGenerator;
    typedef boost::variate_generator<RandomGenerator&, ExponentialDistribution> ExponentialGenerator;
    typedef boost::variate_generator<RandomGenerator&, UniformDistributionInt> UniformGeneratorInt;
    typedef boost::variate_generator<RandomGenerator&, UniformDistributionReal> UniformGeneratorReal;
    typedef boost::variate_generator<RandomGenerator&, PoissonDistribution> PoissonGenerator;

    struct timeval Time;
    gettimeofday(&Time,0);
    static RandomGenerator rng(Time.tv_sec+Time.tv_usec*1000);

    ExponentialDistribution LoiExpoA(1/S);
    UniformDistributionInt LoiUnifInt(0,Npop-1);
    UniformDistributionReal LoiUnifReal(0.,1.);

    ExponentialGenerator TirageExpoA(rng,LoiExpoA);
    UniformGeneratorInt TirageUnifInt(rng,LoiUnifInt);
    UniformGeneratorReal TirageUnifReal(rng,LoiUnifReal);

    int indPop(0);
    for(indPop = 0 ; indPop < Npop ; ++indPop)
    {
        if(FIT == "Derive")
        {
            W.push_back(1.);
        }
        else if(FIT == "Selection")
        {
            Vecteur Selec = SelectionModel1(EChrom1[indPop],SChrom1[indPop],EChrom2[indPop],SChrom2[indPop],h);
            double Wi = Selec.getx1();
            double Hi = Selec.getx2();
            W.push_back(Wi);
            H.push_back(Hi);
            //Test << Wi << endl;
        }
    }

    /*int indMoy(0);
    double MoyW(0);
    for(indMoy = 0 ; indMoy < Npop ; ++indMoy)
    {
        MoyW += W[indMoy];
    }
    MoyW = MoyW/Npop;

    int indPar(0);
    double ParFoc;
    for(indPar = 0 ; indPar < Npop ; ++indPar)
    {
        if(W[indPar] == 0.9)
        {
            ParFoc = indPar;
            break;
        }
    }

    double ComptDel1(0),ComptAa1(0),Comptaa1(0),ComptAA1(0);
    int indDelet(0);
    for(indDelet = 0 ; indDelet < Npop ; ++indDelet)
    {
        if(SChrom1[indDelet] != 0)
        {
            ComptDel1 += 1;
        }
        if(SChrom2[indDelet] != 0)
        {
            ComptDel1 += 1;
        }
        if(SChrom1[indDelet] == 0 && SChrom2[indDelet] == 0)
        {
            ComptAA1 += 1;
        }
        if(SChrom1[indDelet] != 0 && SChrom2[indDelet] == 0)
        {
            ComptAa1 += 1;
        }
        if(SChrom1[indDelet] == 0 && SChrom2[indDelet] != 0)
        {
            ComptAa1 += 1;
        }
        if(SChrom1[indDelet] != 0 && SChrom2[indDelet] != 0)
        {
            Comptaa1 += 1;
        }
    }
    double freqa = ComptDel1/(2*Npop);
    double freqaa = Comptaa1/Npop;
    double freqAa = ComptAa1/Npop;
    double freqAA = ComptAA1/Npop;*/

    //Test << MoyW - (1 - h*S*freqAa - S*freqaa) << endl;

    //Test << ComptDel1/(2*Npop) << endl;

    //vector<double> Repro;

    int ReproSelf,ReproParent1,ReproParent2;
    int indRepro(0);
    double NbreParent(0);
    double Rejetaa(0),RejetAa(0);
    for(indRepro = 0 ; indRepro < Npop ; ++indRepro)
    {
        double ProbaSelfing = TirageUnifReal();
        double ProbaCloning = TirageUnifReal();
        double ProbaAutomixis = TirageUnifReal();

        if(ProbaSelfing < Self)
        {
            //cout << "Self" << endl;
            int Parent;
            do
            {
                Parent = TirageUnifInt();
                double ProbaRepro = TirageUnifReal();
                if(ProbaRepro <= W[Parent])
                {
                    ReproSelf = 1;
                }
                else
                {
                    ReproSelf = 0;
                }
            } while(ReproSelf == 0);

            Vecteur V = RecombinaisonModel1(EChrom1[Parent],SChrom1[Parent],EChrom2[Parent],SChrom2[Parent],EChrom1[Parent],SChrom1[Parent],EChrom2[Parent],SChrom2[Parent],Rjk);
            ECHROM1[indRepro] = V.getx1();
            SCHROM1[indRepro] = V.getx2();
            ECHROM2[indRepro] = V.getx3();
            SCHROM2[indRepro] = V.getx4();
        }
        else if(ProbaCloning < Clone)
        {
            //cout << "Clone" << endl;
            int Parent;
            do
            {
                Parent = TirageUnifInt();
                double ProbaRepro = TirageUnifReal();
                if(ProbaRepro <= W[Parent])
                {
                    ReproSelf = 1;
                }
                else
                {
                    ReproSelf = 0;
                }
            }while(ReproSelf == 0);
            ECHROM1[indRepro] = EChrom1[Parent];
            SCHROM1[indRepro] = SChrom1[Parent];
            ECHROM2[indRepro] = EChrom2[Parent];
            SCHROM2[indRepro] = SChrom2[Parent];
        }
        else if(ProbaAutomixis < Automixis)
        {
            //cout << "Automixis" << endl;
            int Parent;
            do
            {
                //Repro.push_back(1);
                Parent = TirageUnifInt();
                double ProbaRepro = TirageUnifReal();
                //cout << Parent << " " << W[Parent] << endl;
                if(ProbaRepro <= W[Parent])
                {
                    ReproSelf = 1;
                }
                else
                {
                    ReproSelf = 0;
                }
            }while(ReproSelf == 0);
            Vecteur V = RecombinaisonAutomixis(EChrom1[Parent],SChrom1[Parent],EChrom2[Parent],SChrom2[Parent],Rjk,Rij,AUTOMIXIS);
            ECHROM1[indRepro] = V.getx1();
            SCHROM1[indRepro] = V.getx2();
            ECHROM2[indRepro] = V.getx3();
            SCHROM2[indRepro] = V.getx4();
        }
        else
        {
            //Repro.push_back(0);
            //cout << "Normal" << endl;
            int Parent1,Parent2;
            //vector<double> ListRejet1,ListRejet2;

            do
            {
                Parent1 = TirageUnifInt();
                /*double Rejet1;
                do
                {
                    Parent1 = TirageUnifInt();
                    Rejet1 = 0;
                    int i(0);
                    for(i = 0 ; i < ListRejet1.size() ; ++i)
                    {
                        if(ListRejet1[i] == Parent1)
                        {
                            Rejet1 = 1;
                            //cout << "Rejet1" << endl;
                        }
                    }
                } while(Rejet1 == 1);*/

                double ProbaReproParent1 = TirageUnifReal();
                //Test << ProbaReproParent1 << endl;
                /*if(W[Parent1] == 0.9)
                {
                    Test << ProbaReproParent1 << endl;
                        //cout << "1000 !" << endl;
                }*/
                if(ProbaReproParent1 <= W[Parent1])
                {
                    ReproParent1 = 1;
                }
                else
                {
                    ReproParent1 = 0;
                    /*if(W[Parent1] == 0.9)
                    {
                        Rejetaa += 1;
                    }
                    else if (W[Parent1] == 0.975)
                    {
                        RejetAa += 1;
                    }*/
                }
            } while(ReproParent1 == 0);


                    /*do
                    {*/
            do
            {
                Parent2 = TirageUnifInt();
                            /*double Rejet2(0);
                            do
                            {
                                Parent2 = TirageUnifInt();
                                Rejet2 = 0;
                                int i(0);
                                for(i = 0 ; i < ListRejet2.size() ; ++i)
                                {
                                    if(ListRejet2[i] == Parent2)
                                    {
                                        Rejet2 = 1;
                                        //cout << "Rejet 2" <<  endl;
                                    }
                                }
                            } while(Rejet2 == 1);*/

                double ProbaReproParent2 = TirageUnifReal();
                /*if(W[Parent2] == 0.9)
                {
                    Test << ProbaReproParent2 << endl;
                    //cout << "1000 !" << endl;
                }*/

                if(ProbaReproParent2 <= W[Parent2])
                {
                    ReproParent2 = 1;
                }
                else
                {
                    ReproParent2 = 0;
                    /*if(W[Parent2] == 0.9)
                    {
                        Rejetaa += 1;
                    }
                    else if (W[Parent2] == 0.975)
                    {
                        RejetAa += 1;
                    }*/
                                //ListRejet2.push_back(Parent2);
                }
            } while(ReproParent2 == 0);

            /*if(Parent1 == indPar)
            {
                NbreParent = NbreParent + 0.5;
                //cout << "1000 !" << endl;
            }
            if(Parent2 == indPar)
            {
                NbreParent = NbreParent + 0.5;
            }*/


            Vecteur V = RecombinaisonModel1(EChrom1[Parent1],SChrom1[Parent1],EChrom2[Parent1],SChrom2[Parent1],EChrom1[Parent2],SChrom1[Parent2],EChrom2[Parent2],SChrom2[Parent2],Rjk);
            ECHROM1[indRepro] = V.getx1();
            SCHROM1[indRepro] = V.getx2();
            ECHROM2[indRepro] = V.getx3();
            SCHROM2[indRepro] = V.getx4();

            /*double seg1(0),seg2(0),seg3(0),seg4(0);
            seg1=rand() %2;
            seg2=rand() %2;
            if(seg1 == 0)
            {
                ECHROM1[indRepro] = EChrom1[Parent1];
                SCHROM1[indRepro] = SChrom1[Parent1];
            }
            else
            {
                ECHROM1[indRepro] = EChrom2[Parent1];
                SCHROM1[indRepro] = SChrom2[Parent1];
            }
            if(seg2 == 0)
            {
                ECHROM2[indRepro] = EChrom1[Parent2];
                SCHROM2[indRepro] = SChrom1[Parent2];
            }
            else
            {
                ECHROM2[indRepro] = EChrom2[Parent2];
                SCHROM2[indRepro] = SChrom2[Parent2];
            }*/
        }
    }


    //Test << (Rejetaa*(1-freqa))/(RejetAa*freqa) << endl;
    //Test << NbreParent - 0.9/MoyW << endl;
    //cout << NbreParent << endl;

    double ComptDel2(0);
    //double ComptAa2(0),Comptaa2(0),ComptAA2(0);
    int indDelet2(0);
    for(indDelet2 = 0 ; indDelet2 < Npop ; ++indDelet2)
    {
        if(SCHROM1[indDelet2] != 0)
        {
            ComptDel2 += 1;
        }
        if(SCHROM2[indDelet2] != 0)
        {
            ComptDel2 += 1;
        }
        /*if(SCHROM1[indDelet2] == 0 && SCHROM2[indDelet2] == 0)
        {
            ComptAA2 += 1;
        }
        if(SCHROM1[indDelet2] != 0 && SCHROM2[indDelet2] == 0)
        {
            ComptAa2 += 1;
        }
        if(SCHROM1[indDelet2] == 0 && SCHROM2[indDelet2] != 0)
        {
            ComptAa2 += 1;
        }
        if(SCHROM1[indDelet2] != 0 && SCHROM2[indDelet2] != 0)
        {
            Comptaa2 += 1;
        }*/
    }
    /*double freqa2 = ComptDel2/(2*Npop);
    double freqaa2 = Comptaa2/Npop;
    double freqAa2 = ComptAa2/Npop;
    double freqAA2 = ComptAA2/Npop;

    Test << freqa2 - (0.9*freqaa+0.975*0.5*freqAa)/MoyW << endl;*/

    //Test << ComptDel2 - ComptDel1 << endl;
    //Test << ppa - pa*(1-S*pa-0.5*h*S*(1-pa))/MoyW << endl;

    int Na = 2*Npop - ComptDel2;
    //cout << Na << endl;

    if(mutA > 0.)
    {
        PoissonDistribution LoiPoissonA(2*Npop*mutA);
        PoissonGenerator TiragePoissonA(rng,LoiPoissonA);
        int NmutA = TiragePoissonA();
        //Test << NmutA/(2*Npop - ComptDel2) << endl;

        if(NmutA > 0)
        {
            int indMutA(0);
            for(indMutA = 0 ; indMutA < NmutA ; ++indMutA)
            {
                double ChoixChrom;
                int Mutant;

                ChoixChrom = rand() %2 + 1;
                Mutant = TirageUnifInt();
                //double X;

                /*do
                {
                    ChoixChrom = rand() %2 + 1;
                    Mutant = TirageUnifInt();
                    if(ChoixChrom == 1)
                    {
                        X = SCHROM1[Mutant];
                    }
                    else
                    {
                        X = SCHROM2[Mutant];
                    }
                } while(X != 0);*/

                double Mutation;
                if(ALLELES == "2")
                {
                    Mutation = S;
                }
                else if(ALLELES == "Infinite")
                {
                    Mutation = TirageExpoA();
                    //Test << Mutation << endl;
                }

                if(ChoixChrom == 1)
                {
                    SCHROM1[Mutant] = Mutation;
                }
                else if(ChoixChrom == 2)
                {
                    SCHROM2[Mutant] = Mutation;
                }
            }
        }
    }

    /*double ComptDel3(0);
    int indDelet3(0);
    for(indDelet3 = 0 ; indDelet3 < Npop ; ++indDelet3)
    {
        if(SCHROM1[indDelet3] != 0)
        {
            ComptDel3 += 1;
        }
        if(SCHROM2[indDelet3] != 0)
        {
            ComptDel3 += 1;
        }
    }*/

    //Test << ComptDel3/(2*Npop) << endl;

    //Test << (ComptDel3 - ComptDel2)/(2*Npop - ComptDel2) << endl;

    if(mutE > 0.)
    {
        PoissonDistribution LoiPoissonE(2.*Npop*mutE);
        PoissonGenerator TiragePoissonE(rng,LoiPoissonE);
        int NmutE = TiragePoissonE();
        //Test << NmutE << endl;

        if(NmutE > 0)
        {
            int indMutE(0);
            for(indMutE = 0 ; indMutE < NmutE ; ++indMutE)
            {
                double ChoixChrom = rand() %2 + 1;
                int Mutant = TirageUnifInt();
                //cout << Mutant << endl;
                if(ChoixChrom == 1)
                {
                    NormalDistribution LoiNormaleE(0,SigE);
                    GaussianGenerator TirageNormalE(rng,LoiNormaleE);
                    double Mutation = TirageNormalE();
                    //Test << Mutation << endl;
                    ECHROM1[Mutant] = pow(10,log10(ECHROM1[Mutant])+Mutation);

                }
                else if(ChoixChrom == 2)
                {
                    NormalDistribution LoiNormaleE(0,SigE);
                    GaussianGenerator TirageNormalE(rng,LoiNormaleE);
                    double Mutation = TirageNormalE();
                    //Test << Mutation << endl;
                    ECHROM2[Mutant] = pow(10,log10(ECHROM2[Mutant])+Mutation);
                }
            }
        }
    }

    /*double Het(0);
    int indHet(0);
    for(indHet = 0 ; indHet < Npop ; ++indHet)
    {
        if(SCHROM1[indHet] != SCHROM2[indHet])
        {
            Het = Het + 1;
        }
    }
    Het = Het/Npop;*/

    //Test << Het << endl;



    Matrice Descendance(ECHROM1,SCHROM1,ECHROM2,SCHROM2,H,Empty,Empty,Empty,Empty,Empty);

    return Descendance;
}
