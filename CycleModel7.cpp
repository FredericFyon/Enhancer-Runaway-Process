#include "CycleModel7.h"
#include "SelectionModel7.h"
#include "RecombinaisonModel3.h"
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

Matrice CycleModel7(double Ei, Matrice Pop, double mutE1, double mutE2, double mutA, double SigE1, double SigE2, double h, double S, double Rij, double Rjk, double Self, double I, double Concaveness, int Npop, string FIT, string ALLELES, string SHAPE)
{
    vector<double> E1Chrom1 = Pop.getVec1();
    vector<double> E2Chrom1 = Pop.getVec2();
    vector<double> SChrom1 = Pop.getVec3();
    vector<double> E1Chrom2 = Pop.getVec4();
    vector<double> E2Chrom2 = Pop.getVec5();
    vector<double> SChrom2 = Pop.getVec6();
    vector<double> E1CHROM1(Npop),E2CHROM1(Npop),SCHROM1(Npop),E1CHROM2(Npop),E2CHROM2(Npop),SCHROM2(Npop);
    vector<double> W,Empty;

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
            W.push_back(SelectionModel7(Ei,E1Chrom1[indPop],E2Chrom1[indPop],SChrom1[indPop],E1Chrom2[indPop],E2Chrom2[indPop],SChrom2[indPop],h,I));
        }
    }

    if(mutA > 0.)
    {
        PoissonDistribution LoiPoissonA(2.*Npop*mutA);
        PoissonGenerator TiragePoissonA(rng,LoiPoissonA);
        int NmutA = TiragePoissonA();

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

    if(mutE1 > 0.)
    {
        PoissonDistribution LoiPoissonE1(2.*Npop*mutE1);
        PoissonGenerator TiragePoissonE1(rng,LoiPoissonE1);
        int NmutE1 = TiragePoissonE1();

        if(NmutE1 > 0)
        {
            int indMutE1(0);
            for(indMutE1 = 0 ; indMutE1 < NmutE1 ; ++indMutE1)
            {
                double ChoixChrom = rand() %2 + 1;
                int Mutant = TirageUnifInt();
                if(ChoixChrom == 1)
                {
                    NormalDistribution LoiNormaleE1(0,SigE1);
                    GaussianGenerator TirageNormaleE1(rng,LoiNormaleE1);
                    double Mutation = TirageNormaleE1();
                    E1Chrom1[Mutant] = pow(10,log10(E1Chrom1[Mutant])+Mutation);
                }
                else if(ChoixChrom == 2)
                {
                    NormalDistribution LoiNormaleE1(0,SigE1);
                    GaussianGenerator TirageNormaleE1(rng,LoiNormaleE1);
                    double Mutation = TirageNormaleE1();
                    E1Chrom2[Mutant] = pow(10,log10(E1Chrom2[Mutant])+Mutation);
                }
            }
        }
    }

    if(mutE2 > 0.)
    {
        PoissonDistribution LoiPoissonE2(2.*Npop*mutE2);
        PoissonGenerator TiragePoissonE2(rng,LoiPoissonE2);
        int NmutE2 = TiragePoissonE2();

        if(NmutE2 > 0)
        {
            int indMutE2(0);
            for(indMutE2 = 0 ; indMutE2 < NmutE2 ; ++indMutE2)
            {
                double ChoixChrom = rand() %2 + 1;
                int Mutant = TirageUnifInt();
                if(ChoixChrom == 1)
                {
                    NormalDistribution LoiNormaleE2(0,SigE2);
                    GaussianGenerator TirageNormaleE2(rng,LoiNormaleE2);
                    double Mutation = TirageNormaleE2();
                    E2Chrom1[Mutant] = pow(10,log10(E2Chrom1[Mutant])+Mutation);
                }
                else if(ChoixChrom == 2)
                {
                    NormalDistribution LoiNormaleE2(0,SigE2);
                    GaussianGenerator TirageNormaleE2(rng,LoiNormaleE2);
                    double Mutation = TirageNormaleE2();
                    E2Chrom2[Mutant] = pow(10,log10(E2Chrom2[Mutant])+Mutation);
                }
            }
        }
    }

    //cout << "Mutation terminee" << endl;

    int ReproSelf,ReproParent1,ReproParent2;
    int indRepro(0);
    for(indRepro = 0 ; indRepro < Npop ; ++indRepro)
    {
        double ProbaSelfing = TirageUnifReal();
        //cout << ProbaSelfing << endl;
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

            Vecteur V = RecombinaisonModel3(E1Chrom1[Parent],E2Chrom1[Parent],SChrom1[Parent],E1Chrom2[Parent],E2Chrom2[Parent],SChrom2[Parent],E1Chrom1[Parent],E2Chrom1[Parent],SChrom1[Parent],E1Chrom2[Parent],E2Chrom2[Parent],SChrom2[Parent],Rij,Rjk);
            E1CHROM1[indRepro] = V.getx1();
            E2CHROM1[indRepro] = V.getx2();
            SCHROM1[indRepro] = V.getx3();
            E1CHROM2[indRepro] = V.getx4();
            E2CHROM2[indRepro] = V.getx5();
            SCHROM2[indRepro] = V.getx6();
        }
        else
        {
            int Parent1,Parent2;
            do
            {
                Parent1 = TirageUnifInt();
                double ProbaReproParent1 = TirageUnifReal();
                //cout << ProbaReproParent1  << " " << W[Parent1] << endl;
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

            Vecteur V = RecombinaisonModel3(E1Chrom1[Parent1],E2Chrom1[Parent1],SChrom1[Parent1],E1Chrom2[Parent1],E2Chrom2[Parent1],SChrom2[Parent1],E1Chrom1[Parent2],E2Chrom1[Parent2],SChrom1[Parent2],E1Chrom2[Parent2],E2Chrom2[Parent2],SChrom2[Parent2],Rij,Rjk);
            E1CHROM1[indRepro] = V.getx1();
            E2CHROM1[indRepro] = V.getx2();
            SCHROM1[indRepro] = V.getx3();
            E1CHROM2[indRepro] = V.getx4();
            E2CHROM2[indRepro] = V.getx5();
            SCHROM2[indRepro] = V.getx6();
        }
    }

    Matrice Descendance(E1CHROM1,E2CHROM1,SCHROM1,E1CHROM2,E2CHROM2,SCHROM2,Empty,Empty,Empty,Empty);

    return Descendance;
}
