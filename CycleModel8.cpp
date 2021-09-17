#include "CycleModel8.h"
#include "Matrice.h"
#include "Vecteur.h"
#include "SegregationModel8.h"
#include "SelectionModel8.h"
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

Matrice CycleModel8(Matrice Pop, double Ei, double mutT, double mutE, double mutA, double SigT, double SigE, double h, double S, double I, double Rij, double Rjk, double Self, double Clone, int Npop, string FIT, string ALLELES, string SEX)
{
    //ofstream Test("/home/frederic/Bureau/GEFSimulator_CB/bin/Debug/TestMutE1.txt",ios::app);

    vector<double> SexChrom1 = Pop.getVec1();
    vector<double> TChrom1 = Pop.getVec2();
    vector<double> EChrom1 = Pop.getVec3();
    vector<double> SChrom1 = Pop.getVec4();
    vector<double> SexChrom2 = Pop.getVec5();
    vector<double> TChrom2 = Pop.getVec6();
    vector<double> EChrom2 = Pop.getVec7();
    vector<double> SChrom2 = Pop.getVec8();
    vector<double> ECHROM1(Npop),ECHROM2(Npop),SCHROM1(Npop),SCHROM2(Npop),SEXCHROM1(Npop),SEXCHROM2(Npop),TCHROM1(Npop),TCHROM2(Npop);
    vector<double> W,H,Empty;
    double MoyXmale(0),MoyY(0),MoyXfem1(0),MoyXfem2(0);

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

    double NberOfMales(0);
    double NberOfFemales(0);
    int indPop(0);
    for(indPop = 0 ; indPop < Npop ; ++indPop)
    {
        //cout << indPop << endl;
        if(FIT == "Derive")
        {
            W.push_back(1.);
        }
        else if(FIT == "Selection")
        {
            Vecteur Selec = SelectionModel8(Ei,SexChrom1[indPop],TChrom1[indPop],EChrom1[indPop],SChrom1[indPop],SexChrom2[indPop],TChrom2[indPop],EChrom2[indPop],SChrom2[indPop],h,I,SEX);
            double Wi = Selec.getx1();
            //cout << Wi << endl;
            double Hi = Selec.getx2();
            W.push_back(Wi);
            if(Hi != 0)
            {
                 H.push_back(Hi);
            }
            double Xmale = Selec.getx3();
            double Y = Selec.getx4();
            double Xfem1 = Selec.getx5();
            double Xfem2 = Selec.getx6();
            NberOfMales += Selec.getx7();
            NberOfFemales += Selec.getx8();
            //cout << Xmale << " " << Y << " " << Xfem1 << " " << Xfem2 << endl;
            MoyXmale = MoyXmale + Xmale;
            MoyY = MoyY + Y;
            MoyXfem1 = MoyXfem1 + Xfem1;
            MoyXfem2 = MoyXfem2 + Xfem2;
            //Test << Wi << endl;
        }
    }
    //cout << MoyXmale << " " << NberOfMales << endl;
    //cout << MoyXfem1 << " " << NberOfFemales << endl;
    MoyXmale = MoyXmale/NberOfMales;
    MoyY = MoyY/NberOfMales;
    MoyXfem1 = MoyXfem1/NberOfFemales;
    MoyXfem2 = MoyXfem2/NberOfFemales;

    //cout << MoyXmale << " " << MoyY << " " << MoyXfem1 << " " << MoyXfem2 << endl;

    vector<double> MoyExpr(4);
    MoyExpr[0] = MoyXmale;
    MoyExpr[1] = MoyY;
    MoyExpr[2] = MoyXfem1;
    MoyExpr[3] = MoyXfem2;

    //cout << "Fin selection" << endl;

    if(mutA > 0.)
    {
        PoissonDistribution LoiPoissonA(2.*Npop*mutA);
        PoissonGenerator TiragePoissonA(rng,LoiPoissonA);
        int NmutA = TiragePoissonA();
        //Test << NmutA << endl;

        if(NmutA > 0)
        {
            int indMutA(0);
            for(indMutA = 0 ; indMutA < NmutA ; ++indMutA)
            {
                double ChoixChrom = rand() %2 + 1;
                int Mutant = TirageUnifInt();
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
                    SChrom1[Mutant] = Mutation;
                }
                else if(ChoixChrom == 2)
                {
                    SChrom2[Mutant] = Mutation;
                }
            }
        }
    }

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
                if(ChoixChrom == 1)
                {
                    NormalDistribution LoiNormaleE(0,SigE);
                    GaussianGenerator TirageNormalE(rng,LoiNormaleE);
                    double Mutation = TirageNormalE();
                    //cout << Mutation << endl;
                    //Test << Mutation << endl;
                    EChrom1[Mutant] = pow(10,log10(EChrom1[Mutant])+Mutation);

                }
                else if(ChoixChrom == 2)
                {
                    NormalDistribution LoiNormaleE(0,SigE);
                    GaussianGenerator TirageNormalE(rng,LoiNormaleE);
                    double Mutation = TirageNormalE();
                    //cout << Mutation << endl;
                    //Test << Mutation << endl;
                    EChrom2[Mutant] = pow(10,log10(EChrom2[Mutant])+Mutation);
                }
            }
        }
    }

    if(mutT > 0.)
    {
        PoissonDistribution LoiPoissonT(2.*Npop*mutT);
        PoissonGenerator TiragePoissonT(rng,LoiPoissonT);
        int NmutT = TiragePoissonT();
        if(NmutT > 0)
        {
            int indMutT(0);
            for(indMutT = 0 ; indMutT < NmutT ; ++indMutT)
            {
                double ChoixChrom = rand() %2 + 1;
                int Mutant = TirageUnifInt();
                if(ChoixChrom == 1)
                {
                    NormalDistribution LoiNormaleT(0,SigT);
                    GaussianGenerator TirageNormalT(rng,LoiNormaleT);
                    double Mutation = TirageNormalT();
                    //Test << Mutation << endl;
                    TChrom1[Mutant] = pow(10,log10(TChrom1[Mutant])+Mutation);

                }
                else if(ChoixChrom == 2)
                {
                    NormalDistribution LoiNormaleT(0,SigT);
                    GaussianGenerator TirageNormalT(rng,LoiNormaleT);
                    double Mutation = TirageNormalT();
                    //Test << Mutation << endl;
                    TChrom2[Mutant] = pow(10,log10(TChrom2[Mutant])+Mutation);
                }
            }
        }
    }

    //cout << "Fin mutation" << endl;

    int ReproSelf,ReproParent1,ReproParent2;
    int indRepro(0);
    for(indRepro = 0 ; indRepro < Npop ; ++indRepro)
    {
        //cout << TChrom1[indRepro] << " " << EChrom1[indRepro] << " " << SChrom1[indRepro] << " " << TChrom2[indRepro] << " " << EChrom2[indRepro] << " " << SChrom2[indRepro] << endl;
        int Parent1,Parent2;
        do
        {
            Parent1 = TirageUnifInt();
            double ProbaReproParent1 = TirageUnifReal();
            if(ProbaReproParent1 <= W[Parent1])
            {
                ReproParent1 = 1;

                if(SexChrom1[Parent1] == SexChrom2[Parent1])
                {
                     do
                    {
                        do
                        {
                            Parent2 = TirageUnifInt();
                            double ProbaReproParent2 = TirageUnifReal();
                            if((ProbaReproParent2 <= W[Parent2]) && (SexChrom1[Parent2] != SexChrom2[Parent2]))
                            {
                                ReproParent2 = 1;
                            }
                            else
                            {
                                ReproParent2 = 0;
                            }
                        } while(ReproParent2 == 0);
                    } while(Parent2 == Parent1);
                }
                else if(SexChrom1[Parent1] != SexChrom2[Parent1])
                {
                     do
                    {
                        do
                        {
                            Parent2 = TirageUnifInt();
                            double ProbaReproParent2 = TirageUnifReal();
                            if((ProbaReproParent2 <= W[Parent2]) && (SexChrom1[Parent2] == SexChrom2[Parent2]))
                            {
                                ReproParent2 = 1;
                            }
                            else
                            {
                                ReproParent2 = 0;
                            }
                        } while(ReproParent2 == 0);
                    } while(Parent2 == Parent1);
                }
            }
            else
            {
                ReproParent1 = 0;
            }
        } while(ReproParent1 == 0);

        //cout << "Fin choix parents" << endl;

        Vecteur V = SegregationModel8(SexChrom1[Parent1],TChrom1[Parent1],EChrom1[Parent1],SChrom1[Parent1],SexChrom2[Parent1],TChrom2[Parent1],EChrom2[Parent1],SChrom2[Parent1],SexChrom1[Parent2],TChrom1[Parent2],EChrom1[Parent2],SChrom1[Parent2],SexChrom2[Parent2],TChrom2[Parent2],EChrom2[Parent2],SChrom2[Parent2],Rij,Rjk);
        SEXCHROM1[indRepro] = V.getx1();
        TCHROM1[indRepro] = V.getx2();
        ECHROM1[indRepro] = V.getx3();
        SCHROM1[indRepro] = V.getx4();
        SEXCHROM2[indRepro] = V.getx5();
        TCHROM2[indRepro] = V.getx6();
        ECHROM2[indRepro] = V.getx7();
        SCHROM2[indRepro] = V.getx8();
    }

    //cout << "Fin recombinaison" << endl;

    Matrice Descendance(SEXCHROM1,TCHROM1,ECHROM1,SCHROM1,SEXCHROM2,TCHROM2,ECHROM2,SCHROM2,H,MoyExpr);

    return Descendance;
}
