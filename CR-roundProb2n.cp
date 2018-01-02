/*
 *
 *
 *  Created by Maria Abou Chakra on 02/12/14.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include<cstdlib>
#include<fstream>
#include<iostream>
#include<math.h>
#include<map>
#include <string.h>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<vector>
#include<algorithm>
#include<random>


#include "individualCRP-roundProb.h"

using namespace std;

int main(int argc, char* argv[]);
//void outputFunction(std::vector<Individual> &popN, Mymap & phenolist,Mymap & genolist,double *totalC, double *genC,int const& roundN,int const& countGen,int const& inD, int const& oneinD,string optN, int const& gent);
long double factorial(int n);
long double binomialCoeff(int n,int i);
double poissonApprox(int popS, double muta, int k);
int mutateIndividual(std::vector<double> & changeProb, double const & rnNum);
double diffFNs(int const & probOPT, double const& prob, double t, double x,double const& gamma, double const& groupSize);
double weightedFrequency(Individual player, double freqT);
void wrightFisherPopulation(std::vector<Individual> &popN,  int const& NinD, std::mt19937_64 & rng, std::uniform_real_distribution<double> & uniR);
void onegame(std::vector<Individual> &popN0, std::vector<Individual> &popN1,  int const& roundN, int const& inD, int const& gInd,  double const& prob, int const& probOPT,double const& gamma,double *totalC0,double *totalC1,double *totalE0,double *totalE1,double *genC0,double *genC1,string outOpt,  int const & randRound,std::mt19937_64 & rng,std::uniform_real_distribution<double> & uniR,std::uniform_int_distribution<int> & uni2);
void onegeneration (std::vector<Individual> &popN0,std::vector<Individual> &popN1,  int const& roundN, int const& inD, int const& gInd, double const& selectionCo,int const& gameNum,double const& prob, int const& probOPT,double  gamma,double *totalC0,double *totalC1,double *totalE0,double *totalE1,double *genC0,double *genC1, string outOpt, int const & randRound,std::mt19937_64 & rng,std::uniform_real_distribution<double> & uniR, std::uniform_int_distribution<int> & uni2);
void mutatedGeneration(std::vector<Individual> &popN,std::vector<double> & mutationP, int const& inD,int const& roundN,double * sigma,std::mt19937_64 & rng, std::uniform_real_distribution<double> & uniR,  std::uniform_int_distribution<int> & uni2);

int main(int argc, char* argv[])
{
    std::mt19937_64 rng;
    std::random_device sysrand;
    uint32_t random_seed=sysrand();
    rng.seed(random_seed);

 
	int agValue=1;
	//RandNum	RNUm;
	string outOpt(argv[agValue]);agValue++;// is for output opts 1 is for total ave and 2 is a generation file
	int clusterInfo=atoi(argv[agValue]);agValue++;
    int inD=atoi(argv[agValue]);agValue++;
	double selectionCo=atof(argv[agValue]);agValue++;
	int gameNum = atoi(argv[agValue]);agValue++;
	int maxgen=atoi(argv[agValue]);agValue++;///for the looping to a max of genrations
	double mutation=atof(argv[agValue]);agValue++;
	//Game
	int gInd = atoi(argv[agValue]);agValue++;
    double groupSize=double(gInd)+double(gInd);
	int roundN = atoi(argv[agValue]);agValue++;
	double maxPay0=atof(argv[agValue]);agValue++;
    double maxPay1=atof(argv[agValue]);agValue++;
	double indEndowment0=atof(argv[agValue]);agValue++;
    double indEndowment1=atof(argv[agValue]);agValue++;
	double riskEndowment0=atof(argv[agValue]);agValue++;
    double riskEndowment1=atof(argv[agValue]);agValue++;
  	double gameTarget=atof(argv[agValue]);agValue++;
	//Probability functions
	double prob= atof(argv[agValue]);agValue++;
	int probFNOPT=atof(argv[agValue]);agValue++;//this changes it from heavy side, linear to other functions
   	double gamma=atof(argv[agValue]);agValue++;/// is the gamma value in risk probability function
    double sigma[2];
	sigma[0]=atof(argv[agValue]);agValue++; //gaussian mutation, -1 means doesnt exist, for threshold
    sigma[1]=atof(argv[agValue]);agValue++; //gaussian mutation, for contributions
    int randRound=atoi(argv[agValue]);agValue++; //0= all rounds, 1 one random round .... roundN
    double initContribution=atof(argv[agValue]);agValue++; //-1 = random 0-..... is what indiviudals will contribute per round, leads starting condition for the behavior of all individuals
    //int numF=atoi(argv[agValue]);
	double contribution;
	//file names
	ofstream fout;
	
    
    
    std::uniform_int_distribution<int> uniI(0,inD-1);
    std::uniform_real_distribution<double> uniR(0.,1.);
    std::uniform_int_distribution<int> uni2(0,(roundN-1));
    //std::normal_distribution<double> gaus(1,0.05);
	//initializes the population
	//cout<<inD[0]<<" popsize\n";
	
	string nameF ("");
	string loca ("");//location of file
	if(clusterInfo==1)
	{
//
         loca="/home/abouchakra/collective/CRP-roundProb/";
    
	}
    else{
        loca="/Users/abouchakra/Documents/MPI/CollectiveRisk6-CR-new-everyroundprob/CRP-roundProb-program/data/";
    }
	string nameF0 ("");
	string nameF1 ("");
	string nameF2 ("");
    
	//nameF=argv[1];
    
	for(int a=3;a<agValue-1;a++)
	{nameF=nameF+argv[a]+"-";
    //cout<<nameF<<"\n";
    }
	for(int a=3;a<=agValue-1;a++)
		nameF0=nameF0+argv[a]+"-";
	
	char *filN1 = new char[(nameF.length()+205)];
	char *filN2 = new char[(nameF.length()+205)];
	char *filN3a = new char[(nameF.length()+205)];
    char *filN3b = new char[(nameF.length()+205)];
	char *filN4 = new char[(nameF.length()+205)];
	
	double * totalC0 = new double[((roundN)+6)];//caclulates the sum of all the columns over all generations and individuals
    double * totalE0 = new double[((roundN)+1)];//caclulates the sum of all the columns over all generations and individuals
	double * genC0 = new double[(roundN+1)];
    double * totalC1 = new double[((roundN)+6)];//caclulates the sum of all the columns over all generations and individuals
    double * totalE1 = new double[((roundN)+1)];//caclulates the sum of all the columns over all generations and individuals
    double * genC1 = new double[(roundN+1)];
	
	//initializes the population
	std::vector<Individual> popN0(inD);
	std::vector<Individual> popN1(inD);
    
    std::vector<double> mutationProb0(inD+1);
   
    double mp=0;
	int flag=1;
    int mnum=0;
	if(inD>0){
    for(int i2=0;i2<mutationProb0.size();i2++)
	{
		if (flag==1)
		{
			
			if(((inD>=100)&&((double(inD)*mutation)<=10.0))||((inD>=20)&&((mutation)<=0.05)))//condition for Poisson Approx
			{
				mp=mp+poissonApprox(inD,mutation, i2);
			}
			else if (inD<=1000)
			{
				mp=mp+(binomialCoeff(inD,i2)*pow(mutation,(i2))*pow(1-mutation,(inD-i2)));
			}
			else
			{
				std::cout<<"mutation Prob failed\n";
				return 0;//exist main
			}
			if (mp>=1) flag=0;
		}
		if (flag==0){mp=0;}
		mutationProb0[i2]=mp;
	}
    }
 
	//re initialise variables
	int fixed =1;
	int countGen = 0;
	
	for(int i=0;i<((roundN)+6);i++)
	{
		totalC0[i]=0;
        totalC1[i]=0;
        
    }
    for(int i=0;i<((roundN)+1);i++)
    {
        totalE0[i]=0;
        totalE1[i]=0;
        
	}
	//initialize individuals and types in a population
    
    
	for(int i=0;i<inD;i++)
	{
		popN0[i].Initialize(roundN);
		popN0[i].SetTarget(gameTarget);
		popN0[i].SetEndowment(indEndowment0);
        popN0[i].SetOrigEndowment(indEndowment0);
        popN0[i].SetLostEndowment(riskEndowment0);
		popN0[i].SetValues(initContribution, roundN, maxPay0,rng);
		popN0[i].SetIndTarget(groupSize);
        
        popN1[i].Initialize(roundN);
        popN1[i].SetTarget(gameTarget);
        popN1[i].SetEndowment(indEndowment1);
        popN1[i].SetOrigEndowment(indEndowment1);
        popN1[i].SetLostEndowment(riskEndowment1);
        popN1[i].SetValues(initContribution, roundN, maxPay1,rng);
        popN1[i].SetIndTarget(groupSize);
		
	}
 

	
	if (outOpt.find("2",0)!= string::npos)//to output totals of round for each generation
	{
		if(inD>0){
            nameF1= loca+"CRP-0-"+nameF0+"gen.dat";
            strcpy(filN3a, nameF1.c_str());
            fout.open(filN3a);
            fout.close();
            nameF2= loca+"CRP-1-"+nameF0+"gen.dat";
            strcpy(filN3b, nameF2.c_str());
            fout.open(filN3b);
            fout.close();
		}
		
	}
    
    do
    {
        
        for(int i=0;i<(roundN+1);i++)
        {
            if(inD>0)genC0[i]=0;
            if(inD>0)genC1[i]=0;
        }
       // cout<<"before\n";
        onegeneration (popN0,popN1,roundN, inD, gInd, selectionCo,gameNum,prob, probFNOPT, gamma,totalC0,totalC1,totalE0,totalE1,genC0,genC1,outOpt,randRound,rng, uniR,uni2);
       //  cout<<"after\n";
        if (outOpt.find("1",0)!= string::npos)
        {
            
            for(int oneinD=0; oneinD<inD;oneinD++)
            {
                totalC0[0]+=popN0[oneinD].GetOrigEndowment();//average endowment
                totalC0[1]+=popN0[oneinD].GetTarget();//average target
                //cout<<popN[oneinD].GetInfo(2)<<"\n";
                if((popN0[oneinD].GetInfo(2)>0)&& (popN0[oneinD].GetInfo(3)>0)){totalC0[2]+=(popN0[oneinD].GetInfo(3)/(popN0[oneinD].GetInfo(2)*popN0[oneinD].GetOrigEndowment()));}//get average payoff
                if(popN0[oneinD].GetInfo(2)>0){totalC0[3]+=(popN0[oneinD].GetGoalReached()/(popN0[oneinD].GetInfo(2)));}//get goal Reached
                //this collects the frequency of the 4 different strategies C=0, C<R, C=R,C>R
                contribution=popN0[oneinD].CalcContribution();
                totalC0[4]=totalC0[4]+contribution/popN0[oneinD].GetOrigEndowment();
                
                
                totalC1[0]+=popN1[oneinD].GetOrigEndowment();//average endowment
                totalC1[1]+=popN1[oneinD].GetTarget();//average target
                //cout<<popN[oneinD].GetInfo(2)<<"\n";
                if((popN1[oneinD].GetInfo(2)>0)&& (popN1[oneinD].GetInfo(3)>0)){totalC1[2]+=(popN1[oneinD].GetInfo(3)/(popN1[oneinD].GetInfo(2)*popN1[oneinD].GetOrigEndowment()));}//get average payoff
                if(popN1[oneinD].GetInfo(2)>0){totalC1[3]+=(popN1[oneinD].GetGoalReached()/(popN1[oneinD].GetInfo(2)));}//get goal Reached
                //this collects the frequency of the 4 different strategies C=0, C<R, C=R,C>R
                contribution=popN1[oneinD].CalcContribution();
                totalC1[4]=totalC1[4]+contribution/popN1[oneinD].GetOrigEndowment();


            }
            
        }
        
		if (outOpt.find("2",0)!= string::npos)
		{
			
			if(inD>0){
                fout.open(filN3a,ios::app);
                if (!fout)
                {
                    //cout << "Unable to open" << argv[1] << "for appending.\n";
                    fout.open(filN3a); // open for writing
                }
                
                
                for(int i=0;i<(roundN);i++)
                {
                    fout<<genC0[i]/(double(gameNum)*groupSize)<<"  ";
                }
                
                
                fout<<"\n";
                fout.close();
                
                fout.open(filN3b,ios::app);
                if (!fout)
                {
                    //cout << "Unable to open" << argv[1] << "for appending.\n";
                    fout.open(filN3b); // open for writing
                }
                
                
                for(int i=0;i<(roundN);i++)
                {
                    fout<<genC1[i]/(double(gameNum)*groupSize)<<"  ";
                }
                
                
                fout<<"\n";
                fout.close();
			}
           
		}
        
        
        if(inD>0){
            wrightFisherPopulation(popN0, inD, rng,uniR);
            wrightFisherPopulation(popN1, inD, rng,uniR);
            //mutatedGeneration(popN0,mutationProb0, inD,roundN, sigma,RNUm);
            mutatedGeneration(popN0,mutationProb0, inD,roundN, sigma,rng,uniR,uni2);
            mutatedGeneration(popN1,mutationProb0, inD,roundN, sigma,rng,uniR,uni2);
        }
        
        //reset certain values
        
        for(int i=0; i < inD;i++)
        {
            popN0[i].ResetValues();
            popN1[i].ResetValues();
        }
        
        countGen++;
        if(countGen==maxgen)
        {fixed=0;}
        
	}while(fixed == 1);
    
	//cout<<countGen<<"\n";
    if (outOpt.find("1",0)!= string::npos)
    {
        
        if(inD>0)
        {
            nameF1= loca+"CRP-0-"+nameF+"total.dat";
            strcpy(filN1, nameF1.c_str());
            
            fout.open(filN1,ios::app);
            if (!fout)
            {
                //cout << "Unable to open" << argv[1] << "for appending.\n";
                fout.open(filN1); // open for writing
            }
            if(fout.tellp()<=0)
            {
                fout<<"E T payoff targetR contribution ";
                
                for(int id=1; id <=(roundN);id++)
                {
                    fout<<"R"<<id<<" ";
                }
                fout<<"\n";
            }
            
            
            for(int id=0;id<(roundN+5);id++)
            {
                if((id>=5)&&(id<roundN+5))
                {
                    if(totalC0[id]>0){fout<<totalC0[id]/(double(gameNum)*double(gInd)*double(maxgen));}
                    else {fout<<totalC0[id];}
                }
                else{
                    
                    if(totalC0[id]>0){fout<<totalC0[id]/(double(inD)*double(maxgen));}
                    else {fout<<totalC0[id];}
                }
                fout<<"  ";
            }
            fout<<randRound;fout<<"  ";
            fout<<initContribution;
            fout<<"\n";
            
            fout.close();
            
            
            nameF2= loca+"CRP-1-"+nameF+"total.dat";
            strcpy(filN2, nameF2.c_str());
            
            fout.open(filN2,ios::app);
            if (!fout)
            {
                //cout << "Unable to open" << argv[1] << "for appending.\n";
                fout.open(filN2); // open for writing
            }
            if(fout.tellp()<=0)
            {
                fout<<"E T payoff targetR contribution ";
                
                for(int id=1; id <=(roundN);id++)
                {
                    fout<<"R"<<id<<" ";
                }
                fout<<"\n";
            }
            
            
            for(int id=0;id<(roundN+5);id++)
            {
                if((id>=5)&&(id<roundN+5))
                {
                    if(totalC1[id]>0){fout<<totalC1[id]/(double(gameNum)*double(gInd)*double(maxgen));}
                    else {fout<<totalC1[id];}
                }
                else{
                    
                    if(totalC1[id]>0){fout<<totalC1[id]/(double(inD)*double(maxgen));}
                    else {fout<<totalC1[id];}
                }
                fout<<"  ";
            }
            fout<<randRound;fout<<"  ";
            fout<<initContribution;
            fout<<"\n";
            
            fout.close();
            
            
            
        }
    }
    if (outOpt.find("3",0)!= string::npos)
    {
        
        if(inD>0){
            nameF1= loca+"CRP-0-"+nameF+"Endowment.dat";
            strcpy(filN1, nameF1.c_str());
            
            fout.open(filN1,ios::app);
            if (!fout)
            {
                //cout << "Unable to open" << argv[1] << "for appending.\n";
                fout.open(filN1); // open for writing
            }
            if(fout.tellp()<=0)
            {
                
                for(int id=1; id <=(roundN);id++)
                {
                    fout<<"R"<<id<<" ";
                }
                fout<<"\n";
            }
            
            
            for(int id=0;id<(roundN);id++)
            {
                
                if(totalE0[id]>0){  fout<<totalE0[id]/(double(gameNum)*double(gInd)*double(maxgen));}
                else{fout<<totalE0[id];}
                
                fout<<"  ";
            }
            
            fout<<"\n";
            
            fout.close();
            
            nameF2= loca+"CRP-1-"+nameF+"Endowment.dat";
            strcpy(filN2, nameF2.c_str());
            
            fout.open(filN2,ios::app);
            if (!fout)
            {
                //cout << "Unable to open" << argv[1] << "for appending.\n";
                fout.open(filN2); // open for writing
            }
            if(fout.tellp()<=0)
            {
                
                for(int id=1; id <=(roundN);id++)
                {
                    fout<<"R"<<id<<" ";
                }
                fout<<"\n";
            }
            
            
            for(int id=0;id<(roundN);id++)
            {
                
                if(totalE1[id]>0){  fout<<totalE1[id]/(double(gameNum)*double(gInd)*double(maxgen));}
                else{fout<<totalE1[id];}
                
                fout<<"  ";
            }
            
            fout<<"\n";
            
            fout.close();
        }
    
        
    }
    
	///memory release
    delete [] totalC0;
    totalC0 = NULL;
    
    delete [] totalE0;
    totalC0 = NULL;
    
    delete [] genC0;
    genC0 = NULL;
                
                delete [] totalC1;
                totalC0 = NULL;
                
                delete [] totalE1;
                totalC0 = NULL;
                
                delete [] genC1;
                genC0 = NULL;
    
    
	return 0;
}
long double factorial(int n)
{
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
long double binomialCoeff(int n, int i)
{
	return (factorial(n))/(factorial(n-i)*factorial(i));
	
}
double poissonApprox(int popS, double muta, int k)
{
	double lambda=double(popS)*muta;
	return (pow(lambda,k)*exp(-lambda))/factorial(k);
}
//bool memberQ(int *rl, int const& num, int const& listlen)
//{
//	bool opt=false;
//	int m;
//	m = std::count (rl, rl+listlen, num);
//	if(m>0){opt = true;}
//	return opt;
//}
//void getgroup(int *grp, int const& gInd, int const& inD, RandNum	Rnum)
//{
//    
//	int num;
//	for(int i=0; i < gInd;)
//	{
//		num = Rnum.GetRand(0,(inD));
//		grp[i]=num;
//		if (i>0)
//		{
//			if (memberQ(grp, num, (i))== false)//makes sure the same person is not selected twice
//			{i++;}
//			
//		}
//		else {i++;}
//		
//	}
//	return;
//}

double diffFNs(int const & probOPT, double const& prob, double t, double x,double const& gamma, double const& groupSize)
{
     if((probOPT==0)) //Step curve
    {
        return prob;
        
    }
    else if(probOPT==1)
    {
        return (1-pow((x+prob)/(t+prob),gamma)); //Exponential curve 1
    }
    else if(probOPT==2)
    {
        return (prob)*(1.0/((1.0)+exp(gamma*(x-t)))); //Fermi
        
    }
    else if(probOPT==3)
    {
        
        return (prob)*(1.0/((1.0)+exp(gamma*(x-t)/groupSize)));
        
    }
    else if(probOPT==4)
    {
        
        return pow(1.0-(x+prob)/(t+prob),gamma); //Exponential curve 2
        
    }
    else if(probOPT==5)
    {
        return (1.0-((x+prob)/(t+prob)*gamma)); //Piecewise linear risk curve
    }
    else
    {    return 0;
    }
}
double weightedFrequency(Individual player, double freqT)
{
	double wf;
	wf=((player.GetFitness()/freqT));
	return wf;
}
//wright-fisher process, selects based on a frequency or weighted probability
void wrightFisherPopulation(std::vector<Individual> &popN,  int const& NinD, std::mt19937_64 & rng, std::uniform_real_distribution<double> & uniR)
{
	std::vector<Individual> newP(NinD);
	double num;
	double	freqT=0;
	double wfreq[NinD];
	
	for(int i=0; i<NinD; i++)
	{
		freqT = popN[i].GetFitness()+freqT;
	}
	wfreq[0]=weightedFrequency(popN[0], freqT);//start with zero
	for(int i=1; i<NinD; i++)
	{
		wfreq[i]=weightedFrequency(popN[i], freqT)+wfreq[i-1];
	}
	
	for(int n=0; n<NinD;n++)
	{
		//num=RNUm.GetRand01();
		 num=uniR(rng);
		for(int i=0; i<NinD; i++)
		{
			if(num<=(wfreq[i]))
			{
				newP[n]=popN[i];
				i=NinD+1;
			}
		}
	}
	//replace the popN array
	for(int n=0; n<NinD;n++)
	{
		popN[n]=newP[n];
	}
	return;
}
void onegame(std::vector<Individual> &popN0, std::vector<Individual> &popN1,  int const& roundN, int const& inD, int const& gInd,  double const& prob, int const& probOPT,double const& gamma,double *totalC0,double *totalC1,double *totalE0,double *totalE1,double *genC0,double *genC1,string outOpt,  int const & randRound,std::mt19937_64 & rng,std::uniform_real_distribution<double> & uniR,std::uniform_int_distribution<int> & uni2)
{
	//individuals play games
	double num=0;
    double numP[2];
	double totalA=0.0;
	double totalpR=0.0;
	double target;
	double goalY=1.0;
  
//	int group0[gInd];
//    if(inD>0){getgroup(group0, gInd, inD, RNUm);}
//    //for all new games the amount payed goes back to 0

    int oneRandRound=-1;
    
    std::vector<int> myDRvector0;
    std::vector<int> myDRvector1;
    std::vector<int> group0;
     std::vector<int> group1;
    for (int i=0; i<inD; ++i)
    {
        myDRvector0.push_back(i); // 1 2 3 4 5 6 7 8 9 vector of indiviual numbers
        myDRvector1.push_back(i); // 1 2 3 4 5 6 7 8 9 vector of indiviual numbers
        
    }
    std::shuffle (myDRvector0.begin(), myDRvector0.end(),rng); //this shuffles the individuals
    std::shuffle (myDRvector1.begin(), myDRvector1.end(),rng); //this shuffles the individuals
    for (int i=0; i<gInd; ++i)
    {
        group0.push_back(myDRvector0[i]); // 1 2 3 4 5 6 7 8 9 takes the first 5 of the vector of indiviual numbers
       group1.push_back(myDRvector1[i]);
        
    }
    
    
    
    for(int i=0; i < inD;i++)
    {
        
        
        popN0[i].SetInfo(1,0);
        popN1[i].SetInfo(1,0);
        //cout<<i<<" "<<popN0[i].GetEndowment()<<"-"<<popN0[i].GetOrigEndowment()<<"\n";
        popN0[i].SetEndowment(popN0[i].GetOrigEndowment());
        popN1[i].SetEndowment(popN1[i].GetOrigEndowment());
       // cout<<i<<" "<<popN0[i].GetEndowment()<<"-"<<popN0[i].GetOrigEndowment()<<"\n";
    }
  
    
    if(randRound==-2) oneRandRound=uni2(rng); //oneRandRound=RNUm.GetRand(0, roundN); //for each new game this is random but it is the same for one game and only happens once
    
    if(inD>0){target = popN0[0].GetTarget();}
    //else if(inD[1]>0){target = popN0[1].GetTarget();}
    else{target=(gInd);}
    
	for(int r=0; r<roundN; r++)
	{
		//cout<<r<<"\n";
        if(inD>0){
            for(int i=0; i < gInd;i++)
            {
                
                //pop0
                num=popN0[group0[i]].CalcPayment(r,totalpR);//investment simulataneous
                //even if you want to pay something you shouldnt be able to pay more than what you have left
                
                //cout<<group0[i]<<"-i-"; cout<<popN0[group0[i]].GetEndowment()<<"-"<<popN0[group0[i]].GetOrigEndowment()<<"-";cout<<num<<" ";
                
                    if(num >(popN0[group0[i]].GetEndowment()))
                    {num=(popN0[group0[i]].GetEndowment());}
               
              //  cout<<group0[i]<<"-i-"; cout<<popN0[group0[i]].GetEndowment()<<"-"<<popN0[group0[i]].GetOrigEndowment()<<"-";cout<<num<<" ";
                
                popN0[group0[i]].SetBehavior(r, (popN0[group0[i]].GetBehavior(r)+num));//total individual investement per round for all games
              //  cout<<num<<" ";
                if (outOpt.find("1",0)!= string::npos) totalC0[5+r]+=(num);
                   
              //  cout<<totalC0[5+r]<<" ";
                if (outOpt.find("2",0)!= string::npos) genC0[r]+=(num);
             //    cout<<genC0[r]<<" \n";
                popN0[group0[i]].SetInfo(1,(popN0[group0[i]].GetInfo(1)+num));//total individual investment per game
                popN0[group0[i]].SetEndowment(popN0[group0[i]].GetEndowment()-num);
                
            //    cout<<group0[i]<<"-i-"; cout<<popN0[group0[i]].GetEndowment()<<"-"<<popN0[group0[i]].GetOrigEndowment()<<"-";cout<<num<<" ";
             
                totalA = totalA + num;//total payed by the whole group
                
                
                //pop1
                num=popN1[group1[i]].CalcPayment(r,totalpR);//investment simulataneous
                //even if you want to pay something you shouldnt be able to pay more than what you have left
                
                //cout<<group0[i]<<"-i-"; cout<<popN0[group0[i]].GetEndowment()<<"-"<<popN0[group0[i]].GetOrigEndowment()<<"-";cout<<num<<" ";
                
                if(num >(popN1[group1[i]].GetEndowment()))
                {num=(popN1[group1[i]].GetEndowment());}
                
                //  cout<<group0[i]<<"-i-"; cout<<popN0[group0[i]].GetEndowment()<<"-"<<popN0[group0[i]].GetOrigEndowment()<<"-";cout<<num<<" ";
                
                popN1[group1[i]].SetBehavior(r, (popN1[group1[i]].GetBehavior(r)+num));//total individual investement per round for all games
                //  cout<<num<<" ";
                if (outOpt.find("1",0)!= string::npos) totalC1[5+r]+=(num);
                
                //  cout<<totalC0[5+r]<<" ";
                if (outOpt.find("2",0)!= string::npos) genC1[r]+=(num);
                //    cout<<genC0[r]<<" \n";
                popN1[group1[i]].SetInfo(1,(popN1[group1[i]].GetInfo(1)+num));//total individual investment per game
                popN1[group1[i]].SetEndowment(popN1[group1[i]].GetEndowment()-num);
                
                //    cout<<group0[i]<<"-i-"; cout<<popN0[group0[i]].GetEndowment()<<"-"<<popN0[group0[i]].GetOrigEndowment()<<"-";cout<<num<<" ";
                
                totalA = totalA + num;//total payed by the whole group
                
            }
		}
		//	cout<<r<<"\n";
        
        //probability of loosing a part of the endowment based on total group contributions
       // randRound==0 all rounds have a potential risk, randRound==-1 no rounds, randRound==-2 a random round, randRound between 1 and roundN means that one specific round has a potential risk
       
        //individuals types in a population
        
        goalY=1.0;
        
        
        if(totalA<target)//doesnt make sense here but keep it for other versions
        {
            goalY=0;
        }
        
        numP[0]=0;
        numP[1]=0;
        
        if(randRound==-1)//no risk effect right now
        {
            //cout<<"in4\n";
            numP[0]=0;
            numP[1]=0;
        }
        // cout<<randRound<<" "<<r<<" "<<oneRandRound<<"\n";
        if((randRound==0)||((randRound-1)==r)||(oneRandRound==r))
        {
            
            //cout<<randRound<<" "<<r<<" "<<oneRandRound<<"\n";
            
            numP[0]=0;
            numP[1]=0;
//            for (int i=0;i<1;i++)
//            {
//                //cout<<probOPT[i]<<" "<< prob[i]<<" "<< target<<" "<< totalA<<" "<< gamma[i]<<" "<< double(gInd[0])<<"\n";
            
            double probOutput=diffFNs(probOPT, prob, target, totalA, gamma, double(gInd)+double(gInd));
            if (uniR(rng)<=probOutput)
            { numP[0]=1.0;}
           
            
//                if((probOPT==0)&&(totalA>target))
//                {   numP[0]=0;}
//                else if (RNUm.GetRand01()<=diffFNs(probOPT, prob, target, totalA, gamma, double(gInd)))
//                { numP[0]=1.0;}
            
          //  }
        }//this is for controlling which round this happens in
        
        if(inD>0){
            //cout<<"in6\n";
            for(int i=0; i < gInd;i++)
            {
                popN0[group0[i]].SetEndowment(popN0[group0[i]].GetEndowment()*(1.0-popN0[group0[i]].GetLostEndowment()*numP[0]));
                 if (outOpt.find("3",0)!= string::npos) totalE0[r]+=popN0[group0[i]].GetEndowment();
                
                popN1[group1[i]].SetEndowment(popN1[group1[i]].GetEndowment()*(1.0-popN1[group1[i]].GetLostEndowment()*numP[0]));
                if (outOpt.find("3",0)!= string::npos) totalE1[r]+=popN1[group1[i]].GetEndowment();
            }
        }
        
        
		totalpR=totalA;//previous round total
	}
	
    goalY=1.0;
    if(totalA<target)//doesnt make sense here but keep it for other versions
    {
        goalY=0;
    }
	
	//cout<<totalpR<<" "<<target<<" "<<num<<" "<<prob<<"\n";
	//cout<<numP[0]<<" "<<numP[1]<<" "<<prob[0]<<" "<<prob[1]<<"\n";
	if(inD>0){
        
        for(int i=0; i < gInd;i++)
        {
            
            
            //total number of games played
            popN0[group0[i]].SetInfo(2,(popN0[group0[i]].GetInfo(2)+1.0));
            //just adds 1 if goal is met or 0 if not
            //	cout<<goalY<<" "<<popN0[group0[i]].GetGoalReached()<<" ";
            popN0[group0[i]].SetGoalReached(popN0[group0[i]].GetGoalReached()+goalY);
            //	cout<<popN0[group0[i]].GetGoalReached()<<"\n";
            //Payoff calculation of an individual
            popN0[group0[i]].Payoff();
            
            
            
            //total number of games played
            popN1[group1[i]].SetInfo(2,(popN1[group0[i]].GetInfo(2)+1.0));
            //just adds 1 if goal is met or 0 if not
            //	cout<<goalY<<" "<<popN0[group0[i]].GetGoalReached()<<" ";
            popN1[group1[i]].SetGoalReached(popN1[group1[i]].GetGoalReached()+goalY);
            //	cout<<popN0[group0[i]].GetGoalReached()<<"\n";
            //Payoff calculation of an individual
            popN1[group1[i]].Payoff();
            
        }
    }
    return;
}
void onegeneration (std::vector<Individual> &popN0,std::vector<Individual> &popN1,  int const& roundN, int const& inD, int const& gInd, double const& selectionCo,int const& gameNum,double const& prob, int const& probOPT,double  gamma,double *totalC0,double *totalC1,double *totalE0,double *totalE1,double *genC0,double *genC1, string outOpt, int const & randRound,std::mt19937_64 & rng,std::uniform_real_distribution<double> & uniR, std::uniform_int_distribution<int> & uni2)
{
	
	//one generation
	for(int i=0; i < gameNum;i++)
	{
        onegame(popN0,popN1, roundN, inD, gInd, prob, probOPT, gamma,totalC0,totalC1,totalE0,totalE1,genC0,genC1,outOpt,randRound,rng, uniR,uni2);
	}
    
	for(int i=0; i < inD;i++)
	{
		if(popN0[i].GetInfo(2)!=0)
		{
			popN0[i].CalcFitness(selectionCo);
		}
        if(popN1[i].GetInfo(2)!=0)
        {
            popN1[i].CalcFitness(selectionCo);
        }
	}
  
	return;
}
int mutateIndividual(std::vector<double> & changeProb, double const & rnNum)
{
    int grp;
    //	double num;
    //    std::uniform_real_distribution<double> uniR(0.,1.);
    //	num=uniR(rng);
    //std::cout<<n<<" "<<freq.size()<<" i ";
    
    for(int i=0;i<changeProb.size();i++)
    {
        //std::cout<<i<<" ";
        if(rnNum<=changeProb[i])
        {
            grp=i;
            i=changeProb.size()+1;
            break;
        }
    }
    //std::cout<<" g "<<grp<<"\n";
    return grp;
}
void mutatedGeneration(std::vector<Individual> &popN,std::vector<double> & mutationP, int const& inD,int const& roundN,double * sigma,std::mt19937_64 & rng, std::uniform_real_distribution<double> & uniR,  std::uniform_int_distribution<int> & uni2)
{
    int sig=1;
    int maxG=3;
    
    if(roundN==1)maxG=1;
    std::vector<int> myDRvector;
    for (int i=0; i<inD; ++i)
    {
        myDRvector.push_back(i); // 1 2 3 4 5 6 7 8 9 vector of indiviual numbers
       
    }
    std::random_shuffle (myDRvector.begin(), myDRvector.end()); //this shuffles the individuals
    
    int rndIn=0;
    
    for(int i=0; i<maxG; i++)
    {
        int muI=mutateIndividual(mutationP, uniR(rng));
        
        for(int n=0; n<muI;n++)
        {
            //int player=RNUm.GetRand(0,inD);
            int player=myDRvector[rndIn];
			if(i==1){sig=0;}
            else{sig=1;}
            //std::cout<<i<<" "<<sigma[sig]<<"\n";
            //popN[n].SetInfo(0, popN[n].GetInfo(0)+10000);
			popN[player].Mutator(i,uni2(rng),rng, sigma[sig]);
            
            
        }
    }
	return;
}
