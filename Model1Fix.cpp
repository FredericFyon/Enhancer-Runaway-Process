#include <string>
#include <cstring>
#include <vector>
#include <fstream>
#include "Matrice.h"
#include "CycleModel1.h"

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <sys/time.h>

void Model1Fix(double Ei, double SigE, double mutA, double h, double S, double Rjk, double Self, int Npop, int Nfix, int TpsDesinit, string FIT, string ALLELES, string OutputPathway)
{
    typedef boost::uniform_int<int> UniformDistributionInt;
    typedef boost::mt19937 RandomGenerator;
    typedef boost::variate_generator<RandomGenerator&, UniformDistributionInt> UniformGeneratorInt;
    static RandomGenerator rng(clock());

    UniformDistributionInt LoiUnifInt(0,Npop-1);
    UniformGeneratorInt TirUnifInt(rng,LoiUnifInt);

    size_t TailleLineEdit = OutputPathway.size() + 1;
    char* FilePathway = new char[TailleLineEdit];
    strncpy(FilePathway,OutputPathway.c_str(),TailleLineEdit);                                  //Construction du fichier de sortie
    ofstream Results(FilePathway,ios::app);

    vector<double> EChrom1(Npop),SChrom1(Npop),EChrom2(Npop),SChrom2(Npop);                     //Cr�ation des 4 vecteurs repr�sentant les 4 valeurs g�nomiques de chaque individu
    vector<double> Empty;                                                                       //Cr�ation d'un vecteur vide

    int indInitPop(0);
    for(indInitPop = 0 ; indInitPop < Npop ; ++indInitPop)
    {
        EChrom1[indInitPop] = Ei;
        EChrom2[indInitPop] = Ei;
        SChrom1[indInitPop] = 0.;                                                               //Initiation des valeurs g�nomiques : population monomorphique
        SChrom2[indInitPop] = 0.;
    }
    Matrice PopInit(EChrom1,SChrom1,EChrom2,SChrom2,Empty,Empty,Empty,Empty,Empty,Empty);       //Matrice repr�sentant la population initiale                                                 //Cr�ation des matrices stockant les forces des Cis-acting factors et les fitness au cours

    int ComptFix(0);                                                                                            //des g�n�rations et pour les Nit it�rations
    double indFix(0);
    for(indFix = 0 ; indFix < Nfix ; ++indFix)                                                      //Boucle des Nit it�rations de l'�volution d'une m�me population initiale
    {

        Matrice Pop(Empty,Empty,Empty,Empty,Empty,Empty,Empty,Empty,Empty,Empty);               //Cr�ation d'une matrice repr�sentant la population

        int indDesinit(0);
        for(indDesinit = 0 ; indDesinit < TpsDesinit ; ++indDesinit)                            //TpsDesinit g�n�rations sont simul�es pour atteindre l'�quilibre mutation-s�lection
        {                                                                                       //pour la fr�quence en all�le d�l�t�re du g�ne. Pas de mutation au locus du Cis-acting factor
            if(indDesinit == 0)
            {
                Pop = PopInit;                                                                  //Avant que la d�sinitialisation commence, la population est la population initiale
            }
            Pop = CycleModel1(Pop,0.,mutA,SigE,h,S,Rjk,Self,Npop,FIT,ALLELES);
        }

        vector<double> EChrom1_i = Pop.getVec1();
        vector<double> SChrom1_i = Pop.getVec2();
        vector<double> EChrom2_i = Pop.getVec3();
        vector<double> SChrom2_i = Pop.getVec4();

        double chrom = rand() %2 + 1;
        int Mutant = TirUnifInt();
        if(chrom == 1)
        {
            EChrom1_i[Mutant] = Ei + SigE;
        }
        else
        {
            EChrom2_i[Mutant] = Ei + SigE;
        }

        Matrice PopFix(EChrom1_i,SChrom1_i,EChrom2_i,SChrom2_i,Empty,Empty,Empty,Empty,Empty,Empty);

        int NbAllMut = 0;

        do
        {
            PopFix = CycleModel1(PopFix,0.,mutA,SigE,h,S,Rjk,Self,Npop,FIT,ALLELES);

            vector<double> EChrom1_Fix = Pop.getVec1();
            vector<double> EChrom2_Fix = Pop.getVec3();

            int indPop(0);
            for(indPop = 0 ; indPop < Npop ; ++indPop)
            {
                if((EChrom1_Fix[indPop] != Ei) && (EChrom2_Fix[indPop] == Ei))
                {
                    NbAllMut += 1;
                }
                else if((EChrom1_Fix[indPop] == Ei) && (EChrom2_Fix[indPop] != Ei))
                {
                    NbAllMut += 1;
                }
                else if((EChrom1_Fix[indPop] != Ei) && (EChrom2_Fix[indPop] != Ei))
                {
                    NbAllMut += 2;
                }
            }
        } while((NbAllMut != 0) && (NbAllMut != 2*Npop));

        if(NbAllMut == 2*Npop)
        {
            ComptFix += 1;
        }
    }

    if(Results)
    {
        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        Results << "RESULTATS" << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        Results << "Valeurs des param�tres : " << endl;
        Results << endl;
        Results << "Force initiale des Cis-acting factors : " << Ei << endl;
        Results << "Force des Cis-acting factors mutants" << Ei+SigE << endl;
        Results << "Coefficient de dominance : " << h << endl;
        Results << "Intensit� de la s�lection au Locus A : " << S << endl;
        Results << "Taux de mutation au Locus A : " << mutA << endl;
        Results << "Taux de recombinaison entre les deux locus : " << Rjk << endl;
        Results << "Taux d'autof�condation : " << Self << endl;
        Results << "Taille de la population d'individus diplo�des : " << Npop << endl;
        Results << "Nombre d'it�rations de simulations de fixation : " << Nfix << endl;
        Results << "Nombre de g�n�rations simul�es pour d�sinitialiser la population : " << TpsDesinit << endl;
        if(FIT == "Selection")
        {
            Results << "S�lection active" << endl;
        }
        else if(FIT == "Derive")
        {
            Results << "D�rive" << endl;
        }
        if(ALLELES == "2")
        {
            Results << "2 all�les au Locus A" << endl;
        }
        else if(ALLELES == "Infinite")
        {
            Results << "Nombre infini d'all�les au locus A" << endl;
        }
        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        Results << "Nombre de fixations de l'all�le mutant observ�es : " << ComptFix << endl;
    }
}
