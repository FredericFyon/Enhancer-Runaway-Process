#include "CycleModel4.h"
#include "SelectionModel4.h"
#include "RecombinaisonModel4.h"
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

Matrice CycleModel4(Matrice Pop, double Ei, double mutR, double mutE, double mutA, double SigR, double SigE, double h, double S, double Rij, double Self, double I, int Npop, string FIT, string ALLELES)
{
    vector<double> RChrom1 = Pop.getVec1();
    vector<double> EChrom1 = Pop.getVec2();
    vector<double> SChrom1 = Pop.getVec3();
    vector<double> RChrom2 = Pop.getVec4();
    vector<double> EChrom2 = Pop.getVec5();
    vector<double> SChrom2 = Pop.getVec6();
    vector<double> RCHROM1(Npop),ECHROM1(Npop),SCHROM1(Npop),RCHROM2(Npop),ECHROM2(Npop),SCHROM2(Npop);
    vector<double> W,Empty,H;

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
            Vecteur Selec = SelectionModel4(Ei,EChrom1[indPop],SChrom1[indPop],EChrom2[indPop],SChrom2[indPop],h,I);
            W.push_back(Selec.getx1());
            double hi = Selec.getx2();
            if(hi != 0)
            {
                H.push_back(hi);
            }
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
                else
                {
                    SChrom2[Mutant] = Mutation;
                }
            }
        }
    }

    //cout << mutR << endl;

    if(mutR > 0.)
    {
        PoissonDistribution LoiPoissonR(2.*Npop*mutR);
        PoissonGenerator TiragePoissonR(rng,LoiPoissonR);
        int NmutR = TiragePoissonR();
        //cout << NmutR << endl;

        if(NmutR > 0)
        {
            int indMutR(0);
            for(indMutR = 0 ; indMutR < NmutR ; ++indMutR)
            {
                double ChoixChrom = rand() %2 + 1;
                int Mutant = TirageUnifInt();
                if(ChoixChrom == 1)
                {
                    NormalDistribution LoiNormaleR(0,SigR);
                    GaussianGenerator TirageNormaleR(rng,LoiNormaleR);
                    double Mutation = TirageNormaleR();
                    RChrom1[Mutant] = pow(10,log10(RChrom1[Mutant])+Mutation);
                }
                else if(ChoixChrom == 2)
                {
                    NormalDistribution LoiNormaleR(0,SigR);
                    GaussianGenerator TirageNormaleR(rng,LoiNormaleR);
                    double Mutation = TirageNormaleR();
                    RChrom2[Mutant] = pow(10,log10(RChrom2[Mutant])+Mutation);
                }
            }
        }
    }

    if(mutE > 0.)
    {
        PoissonDistribution LoiPoissonE(2.*Npop*mutE);
        PoissonGenerator TiragePoissonE(rng,LoiPoissonE);
        int NmutE = TiragePoissonE();

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
                    GaussianGenerator TirageNormaleE(rng,LoiNormaleE);
                    double Mutation = TirageNormaleE();
                    EChrom1[Mutant] = pow(10,log10(EChrom1[Mutant])+Mutation);
                }
                else if(ChoixChrom == 2)
                {
                    NormalDistribution LoiNormaleE(0,SigE);
                    GaussianGenerator TirageNormaleE(rng,LoiNormaleE);
                    double Mutation = TirageNormaleE();
                    EChrom2[Mutant] = pow(10,log10(EChrom2[Mutant])+Mutation);
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

            Vecteur V = RecombinaisonModel4(RChrom1[Parent],EChrom1[Parent],SChrom1[Parent],RChrom2[Parent],EChrom2[Parent],SChrom2[Parent],RChrom1[Parent],EChrom1[Parent],SChrom1[Parent],RChrom2[Parent],EChrom2[Parent],SChrom2[Parent],Rij);
            RCHROM1[indRepro] = V.getx1();
            ECHROM1[indRepro] = V.getx2();
            SCHROM1[indRepro] = V.getx3();
            RCHROM2[indRepro] = V.getx4();
            ECHROM2[indRepro] = V.getx5();
            SCHROM2[indRepro] = V.getx6();
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

            Vecteur V = RecombinaisonModel4(RChrom1[Parent1],EChrom1[Parent1],SChrom1[Parent1],RChrom2[Parent1],EChrom2[Parent1],SChrom2[Parent1],RChrom1[Parent2],EChrom1[Parent2],SChrom1[Parent2],RChrom2[Parent2],EChrom2[Parent2],SChrom2[Parent2],Rij);
            RCHROM1[indRepro] = V.getx1();
            ECHROM1[indRepro] = V.getx2();
            SCHROM1[indRepro] = V.getx3();
            RCHROM2[indRepro] = V.getx4();
            ECHROM2[indRepro] = V.getx5();
            SCHROM2[indRepro] = V.getx6();
        }
    }

    Matrice Descendance(RCHROM1,ECHROM1,SCHROM1,RCHROM2,ECHROM2,SCHROM2,H,Empty,Empty,Empty);

    return Descendance;
}
