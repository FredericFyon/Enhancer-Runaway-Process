#include "CycleModel2.h"
#include "SelectionModel2.h"
#include "RecombinaisonModel2.h"
#include "Matrice.h"
#include "Vecteur.h"
#include <string>

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
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>

using namespace std;

Matrice CycleModel2(Matrice Pop, double mutE, double mutA, double SigE, double h, double S, double Rjk, double Self, double I, double Concaveness, int Npop, string FIT, string ALLELES, string SHAPE)
{
    vector<double> EChrom1a = Pop.getVec1();
    vector<double> SChrom1a = Pop.getVec2();
    vector<double> EChrom2a = Pop.getVec3();
    vector<double> SChrom2a = Pop.getVec4();
    vector<double> EChrom1b = Pop.getVec5();
    vector<double> SChrom1b = Pop.getVec6();
    vector<double> EChrom2b = Pop.getVec7();
    vector<double> SChrom2b = Pop.getVec8();
    vector<double> W,Empty;
    vector<double> ECHROM1A(Npop),SCHROM1A(Npop),ECHROM2A(Npop),SCHROM2A(Npop),ECHROM1B(Npop),SCHROM1B(Npop),ECHROM2B(Npop),SCHROM2B(Npop);

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
            W.push_back(SelectionModel2(EChrom1a[indPop],SChrom1a[indPop],EChrom2a[indPop],SChrom2a[indPop],EChrom1b[indPop],SChrom1b[indPop],EChrom2b[indPop],SChrom2b[indPop],h,I,Concaveness,SHAPE));
        }
    }

    if(mutA > 0.)
    {
        PoissonDistribution LoiPoissonA(2.*Npop*mutA);
        PoissonGenerator TiragePoissonA(rng,LoiPoissonA);
        int NmutAa = TiragePoissonA();
        int NmutAb = TiragePoissonA();

        if(NmutAa > 0)
        {
            int indMutAa(0);
            for(indMutAa = 0 ; indMutAa < NmutAa ; ++indMutAa)
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
                }

                if(ChoixChrom == 1)
                {
                    SChrom1a[Mutant] = Mutation;
                }
                else if(ChoixChrom == 2)
                {
                    SChrom2a[Mutant] = Mutation;
                }
            }
        }

        if(NmutAb > 0)
        {
            int indMutAb(0);
            for(indMutAb = 0 ; indMutAb < NmutAb ; ++indMutAb)
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
                }

                if(ChoixChrom == 1)
                {
                    SChrom1b[Mutant] = Mutation;
                }
                else if(ChoixChrom == 2)
                {
                    SChrom2b[Mutant] = Mutation;
                }
            }
        }
    }

    if(mutE > 0.)
    {
        PoissonDistribution LoiPoissonE(2.*Npop*mutE);
        PoissonGenerator TiragePoissonE(rng,LoiPoissonE);
        int NmutEa = TiragePoissonE();
        int NmutEb = TiragePoissonE();

        if(NmutEa > 0)
        {
            int indMutE(0);
            for(indMutE = 0 ; indMutE < NmutEa ; ++indMutE)
            {
                double ChoixChrom = rand() %2 + 1;
                int Mutant = TirageUnifInt();
                if(ChoixChrom == 1)
                {
                    NormalDistribution LoiNormaleE(0,SigE);
                    GaussianGenerator TirageNormalE(rng,LoiNormaleE);
                    double Mutation = TirageNormalE();
                    EChrom1a[Mutant] = pow(10,log10(EChrom1a[Mutant])+Mutation);
                }
                else if(ChoixChrom == 2)
                {
                    NormalDistribution LoiNormaleE(0,SigE);
                    GaussianGenerator TirageNormalE(rng,LoiNormaleE);
                    double Mutation = TirageNormalE();
                    EChrom2a[Mutant] = pow(10,log10(EChrom2a[Mutant])+Mutation);
                }
            }
        }

        if(NmutEb > 0)
        {
            int indMutE(0);
            for(indMutE = 0 ; indMutE < NmutEb ; ++indMutE)
            {
                double ChoixChrom = rand() %2 + 1;
                int Mutant = TirageUnifInt();
                if(ChoixChrom == 1)
                {
                    NormalDistribution LoiNormaleE(0,SigE);
                    GaussianGenerator TirageNormalE(rng,LoiNormaleE);
                    double Mutation = TirageNormalE();
                    EChrom1b[Mutant] = pow(10,log10(EChrom1b[Mutant])+Mutation);
                }
                else if(ChoixChrom == 2)
                {
                    NormalDistribution LoiNormaleE(0,SigE);
                    GaussianGenerator TirageNormalE(rng,LoiNormaleE);
                    double Mutation = TirageNormalE();
                    EChrom2b[Mutant] = pow(10,log10(EChrom2b[Mutant])+Mutation);
                }
            }
        }
    }

    int ReproSelf,ReproParent1,ReproParent2;
    int indRepro(0);
    for(indRepro = 0 ; indRepro < Npop ; ++indRepro)
    {
        double ProbaSelfing = TirageUnifReal();
        if(ProbaSelfing <= Self)
        {
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

            Vecteur V = RecombinaisonModel2(EChrom1a[Parent],SChrom1a[Parent],EChrom2a[Parent],SChrom2a[Parent],EChrom1b[Parent],SChrom1b[Parent],EChrom2b[Parent],SChrom2b[Parent],EChrom1a[Parent],SChrom1a[Parent],EChrom2a[Parent],SChrom2a[Parent],EChrom1b[Parent],SChrom1b[Parent],EChrom2b[Parent],SChrom2b[Parent],Rjk);
            ECHROM1A[indRepro] = V.getx1();
            SCHROM1A[indRepro] = V.getx2();
            ECHROM2A[indRepro] = V.getx3();
            SCHROM2A[indRepro] = V.getx4();
            ECHROM1B[indRepro] = V.getx5();
            SCHROM1B[indRepro] = V.getx6();
            ECHROM2B[indRepro] = V.getx7();
            SCHROM2B[indRepro] = V.getx8();
        }
        else
        {
            int Parent1,Parent2;
            do
            {
                Parent1 = TirageUnifInt();
                double ProbaReproParent1 = TirageUnifReal();
                if(ProbaReproParent1 <= W[Parent1])
                {
                    ReproParent1 = 1;

                    do
                    {
                        do
                        {
                            Parent2 = TirageUnifInt();
                            double ProbaReproParent2 = TirageUnifReal();
                            if(ProbaReproParent2 <= W[Parent2])
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
                else
                {
                    ReproParent1 = 0;
                }
            } while(ReproParent1 == 0);

            Vecteur V = RecombinaisonModel2(EChrom1a[Parent1],SChrom1a[Parent1],EChrom2a[Parent1],SChrom2a[Parent1],EChrom1b[Parent1],SChrom1b[Parent1],EChrom2b[Parent1],SChrom2b[Parent1],EChrom1a[Parent2],SChrom1a[Parent2],EChrom2a[Parent2],SChrom2a[Parent2],EChrom1b[Parent2],SChrom1b[Parent2],EChrom2b[Parent2],SChrom2b[Parent2],Rjk);
            ECHROM1A[indRepro] = V.getx1();
            SCHROM1A[indRepro] = V.getx2();
            ECHROM2A[indRepro] = V.getx3();
            SCHROM2A[indRepro] = V.getx4();
            ECHROM1B[indRepro] = V.getx5();
            SCHROM1B[indRepro] = V.getx6();
            ECHROM2B[indRepro] = V.getx7();
            SCHROM2B[indRepro] = V.getx8();
        }
    }

    Matrice Descendance(ECHROM1A,SCHROM1A,ECHROM2A,SCHROM2A,ECHROM1B,SCHROM1B,ECHROM2B,SCHROM2B,Empty,Empty);

    return Descendance;
}
