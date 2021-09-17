#include <fstream>
#include <iostream>
#include <string>
#include "Model1.h"
#include "Model1Fix.h"
#include "Model2.h"
#include "Model3.h"
#include "Model4.h"
#include "Model5.h"
#include "Model6.h"
#include "Model7.h"
#include "Model8.h"
#include "ModelRFix.h"
#include "ModelDom2All.h"

using namespace std;

int main()
{
    ifstream param("/home/ffyon/GEFSimulator/Param.txt");
    //ifstream param("/home/frederic/Bureau/GEFSimulator_CB/bin/Debug/Param.txt");
    double Ei, mutT,mutE,mutA,SigT,SigE,h,S,Rij,Rjk,Self,Clone,Automixis,Concaveness,I;
    int Npop,Ngen,Ngen2,Nit,TpsDesinit,PasEchant;
    string ModelType,FIT,ALLELES,SUIVIFIT,OutputPathway,SHAPE,Log,CANCER,HOMOZYGOSITE,AUTOMIXIS,SEX;

    param.seekg(191,ios::cur);
    param >> ModelType;
    //cout << ModelType << endl;
    param >> Ei;
    //cout << Ei << endl;
    param >> mutT;
    //cout << mutT << endl;
    param >> mutE;
    //cout << mutE << endl;
    param >> mutA;
    //cout << mutA << endl;
    param >> SigT;
    //cout << SigT << endl;
    param >> SigE;
    //cout << SigE << endl;
    param >> h;
    //cout << h << endl;
    param >> S;
    //cout << S << endl;
    param >> Rij;
    //cout << Rij << endl;
    param >> Rjk;
    //cout << Rjk << endl;
    param >> Self;
    param >> Clone;
    param >> Automixis;
    param >> Concaveness;
    param >> I;
    param >> Npop;
    param >> Ngen;
    param >> Ngen2;
    param >> Nit;
    param >> TpsDesinit;
    param >> PasEchant;
    param >> FIT;
    param >> ALLELES;
    param >> SUIVIFIT;
    param >> SHAPE;
    param >> Log;
    param >> CANCER;
    param >> HOMOZYGOSITE;
    param >> AUTOMIXIS;
    param >> SEX;
    param >> OutputPathway;
    //cout << OutputPathway << endl;

    if(ModelType == "BasicModel")
    {
        Model1(Ei,mutE,mutA,SigE,h,S,Rij,Rjk,Self,Clone,Automixis,Npop,Ngen,Nit,TpsDesinit,PasEchant,FIT,ALLELES,SUIVIFIT,OutputPathway,Log,CANCER,HOMOZYGOSITE,AUTOMIXIS);
    }
    /*else if(ModelType == "BasicFixModel")
    {
        Model1Fix(Ei,SigE,mutA,h,S,Rjk,Self,Npop,Nit,TpsDesinit,FIT,ALLELES,OutputPathway);
    }*/
    else if(ModelType == "DuplicatedModel")
    {
        Model2(Ei,mutE,mutA,SigE,h,S,Rjk,Self,Concaveness,I,Npop,Ngen,Nit,TpsDesinit,PasEchant,FIT,ALLELES,SUIVIFIT,OutputPathway,SHAPE);
    }
    else if(ModelType == "CisTransModel")
    {
        Model3(Ei,mutT,mutE,mutA,SigT,SigE,h,S,Rij,Rjk,Self,I,Concaveness,Npop,Ngen,Nit,TpsDesinit,PasEchant,FIT,ALLELES,SUIVIFIT,OutputPathway,SHAPE);
    }
    else if(ModelType == "RModifierModel")
    {
        Model4(Ei,mutT,mutE,mutA,SigT,SigE,h,S,Rij,Rjk,Self,I,Npop,Ngen,Ngen2,Nit,TpsDesinit,PasEchant,FIT,ALLELES,SUIVIFIT,OutputPathway,HOMOZYGOSITE);
    }
    else if(ModelType == "DModifierModel")
    {
        Model5(Ei,mutT,mutE,mutA,SigT,SigE,h,S,Rij,Rjk,Self,I,Npop,Ngen,Ngen2,Nit,TpsDesinit,PasEchant,FIT,ALLELES,SUIVIFIT,OutputPathway,HOMOZYGOSITE);
    }
    else if(ModelType == "QProtModifierModel")
    {
        Model6(Ei,mutT,mutE,mutA,SigT,SigE,h,S,Rij,Rjk,Self,I,Npop,Ngen,Ngen2,Nit,TpsDesinit,PasEchant,FIT,ALLELES,SUIVIFIT,OutputPathway,"QProt",Log);
    }
    else if(ModelType == "ExcesQProtModifierModel")
    {
        Model6(Ei,mutT,mutE,mutA,SigT,SigE,h,S,Rij,Rjk,Self,I,Npop,Ngen,Ngen2,Nit,TpsDesinit,PasEchant,FIT,ALLELES,SUIVIFIT,OutputPathway,"ExcesQProt",Log);
    }
    else if(ModelType == "CisCisModifierModel")
    {
        Model7(Ei,mutT,mutE,mutA,SigT,SigE,h,S,Rij,Rjk,Self,I,Npop,Ngen,Ngen2,Nit,TpsDesinit,PasEchant,Concaveness,FIT,ALLELES,SUIVIFIT,OutputPathway,SHAPE);
    }
    else if(ModelType == "SexModel")
    {
        Model8(Ei,mutT,mutE,mutA,SigT,SigE,h,S,I,Rij,Rjk,Self,Clone,Npop,Ngen,Nit,TpsDesinit,PasEchant,FIT,ALLELES,SUIVIFIT,OutputPathway,Log,SEX);
    }
    /*else if(ModelType == "FixationR")
    {
        ModelRFix(Ei,mutE,mutA,SigT,SigE,h,S,Rij,Rjk,Self,I,Npop,Nit,TpsDesinit,FIT,ALLELES,SUIVIFIT,OutputPathway);
    }*/
    else if(ModelType == "DModifier2AllModel")
    {
        ModelDom2All(Ei,mutT,mutE,mutA,SigT,SigE,h,S,Rij,Rjk,Self,I,Npop,Ngen,Ngen2,Nit,TpsDesinit,PasEchant,FIT,ALLELES,SUIVIFIT,OutputPathway);
    }
    else
    {
        cout << "Model type has not been well defined. Please try again using a valid reference (<<BasicModel>>, <<DuplicatedModel>>, <<CisTransModel>> or <<RModifierModel>>)" << endl;
    }
}
