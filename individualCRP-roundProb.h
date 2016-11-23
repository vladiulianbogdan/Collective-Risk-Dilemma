/*
 *  individual.h
 *  
 *
 *  Created by Maria Abou Chakra on 2/28/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

//#include<iostream>
//using namespace std;
#include<vector>
class Individual
{
public:
	Individual();//constructor fn
	void Initialize(int roundN);//constructor fn
	void ReInitialize(int roundN);//constructor fn
	~Individual();//destructor fn
	void SetValues(double initnum, int roundN, double maxP,std::mt19937_64 & rng);
	void ResetValues();
	double GetInfo(int r); //accessor fn
	double SetInfo(int r, double value); //accessor fn
	double GetBehavior(int r);
	double SetBehavior(int r,double value);
	double Getgenebelow(int r);
	double Getgeneabove(int r);
	double GetT1(int r);
	double Payoff(double &prob);
    double Payoff();
	void CalcFitness(double const& selectionCo);
	void SetFitness(double value);
	double GetFitness();
	void SetGoalReached(double value);
	double GetGoalReached();
	int GetArraySize();
	void SetEndowment(double value);
	double CalcContribution();
	double GetEndowment();
    void SetOrigEndowment(double value);
    double GetOrigEndowment();
    void SetLostEndowment(double value);
    double GetLostEndowment();
	void SetTarget(double value);
	void SetIndTarget(double value);
	double GetTarget();
	double GetIndTarget();
	void SetMaxPayment(double value);
	double GetMaxPayment();
	void Mutator(int mp, int r,std::mt19937_64 & rng, double sigma);
	double CalcPayment(int roundN, double const & totalPotsize);
private:
	double ind[4];
	std::vector<double> behaviorR;
	std::vector<double> genethreshold1;
	std::vector<double> genebelow;
	std::vector<double> geneabove;
	double indFit;
	double endowment;
    double originalEndowment;
	double contribution;
	double maxpayment;
	double target;
	double individualtarget;
    double endowmentLost;
	double indGoalR;
	int info;
	int acol;
	
};

//constructor
Individual::Individual()
{
	info=4;
	//ind= new int*[info];
	/*
	 ind[0] = Number
	 ind[1] = onegame
	 ind[2] = games played
	 ind[3] = payoff
	 behaviorR = round payments one for each round
	 genebelow/geneabove/genemiddle=how much one would pay/thresshold
	 genethreshold1/genethreshold2= threshold*/
	indFit=1.0;
	endowment=0;
	target=0;
	indGoalR=0;
	
}
//constructor
void Individual::Initialize(int roundN)
{
	info=4;//arraysize for individual information
	acol = roundN;
	behaviorR.reserve(roundN);
	genethreshold1.reserve(roundN);
	genebelow.reserve(roundN);
	geneabove.reserve(roundN);
}
//constructor
void Individual::ReInitialize(int roundN)
{
	info=4;//arraysize for individual information
	acol = roundN;
	behaviorR.clear();
	genethreshold1.clear();
	genebelow.clear();
	geneabove.clear();
	
}
//Destructor
Individual::~Individual()
{

}
void Individual::SetValues(double initnum, int roundN,double maxP, std::mt19937_64 & rng)
{
	maxpayment=maxP;
	indFit=1.0;
	
    std::uniform_real_distribution<double> uniR(0.,1.);
		
	ind[0]=initnum;//individual Number
	ind[1]=0;//amount payed each game
	ind[2]=0;//total number of games played
	ind[3]=0;//payoff total
	for(int i=0; i < (roundN);i++)
	{ 
		
		behaviorR.push_back(0);//initialize the action for each round
        if(i==0)
        {
            genethreshold1.push_back(0);
            if (initnum<0)
            {
                genebelow.push_back(uniR(rng)*maxpayment);
            }
            else{
                genebelow.push_back(initnum);
            }//genebelow.push_back(maxpayment);
            //genebelow.push_back(0);
            
            geneabove.push_back(genebelow[i]);
        }
        else{
            genethreshold1.push_back(uniR(rng));
            if (initnum<0)
            {
                genebelow.push_back(uniR(rng)*maxpayment);
                geneabove.push_back(uniR(rng)*maxpayment);
            }
            else{
                genebelow.push_back(initnum);
                geneabove.push_back(initnum);
            }
        }
	/*std::cout<<genethreshold1[i]<<" ";
		std::cout<<genebelow[i]<<" ";
		std::cout<<geneabove[i]<<" ";*/
	
	 }
	//std::cout<<"\n";
		return;
}
void Individual::ResetValues()
{
	
	indFit=1.0;
	indGoalR=0;
	ind[1]=0;//amount payed each game
	ind[2]=0;//total number of games played
	ind[3]=0;//payoff total
	
	for(int i=0; i < (acol);i++)
	{ 	
		behaviorR[i]=0;//initialize the action for each round
	}
	return;
}
int Individual::GetArraySize()
{
	return info;
}
double Individual::GetInfo(int r)
{
	return ind[r];
}
double Individual::SetInfo(int r, double value)
{
	ind[r] = value;
	return ind[r];
}
double Individual::GetBehavior(int r)
{
	return behaviorR[r];
}
double Individual::SetBehavior(int r,double value)
{
	behaviorR[r] = value;
	return behaviorR[r];
}
double Individual::Getgenebelow(int r)
{
	return genebelow[r];
}
double Individual::Getgeneabove(int r)
{
	return geneabove[r];
}
double Individual::GetT1(int r)
{
	return genethreshold1[r];
}
double Individual::Payoff(double &prob)
{
	ind[3]=ind[3]+((endowment-ind[1])*prob);
	return ind[3];
}
double Individual::Payoff()
{
    ind[3]=ind[3]+(endowment);
    return ind[3];
}
void Individual::CalcFitness(double const& selectionCo)
{
	
	indFit=exp(selectionCo*(double(ind[3])/double(ind[2])));
	return;
}
void Individual::SetFitness(double value)
{
	indFit=value;
	return;
}
double Individual::GetFitness()
{
	
	return indFit;
}
void Individual::SetGoalReached(double value)
{
	indGoalR = value;
	return;
}
double Individual::GetGoalReached()
{
	return indGoalR;
}
void Individual::SetEndowment(double value)
{
	endowment=value;
	return;
}
void Individual::SetOrigEndowment(double value)
{
    originalEndowment=value;
    return;
}
double Individual::GetOrigEndowment()
{
    return originalEndowment;
}
void Individual::SetLostEndowment(double value)
{
    endowmentLost=value;
    return;
}
double Individual::GetLostEndowment()
{
    return endowmentLost;
}
double Individual::CalcContribution()
{
	contribution=0;
	for(int r=0;r<acol;r++)
		contribution+=(behaviorR[r]);
	
	if(ind[2]>0)return contribution/(ind[2]);
	else{return 1.0;}
}
double Individual::GetEndowment()
{
	return endowment;
}
void Individual::SetTarget(double value)
{
	target=value;
	return;
}
void Individual::SetIndTarget(double value)
{
	individualtarget=target/value;
	return;
}
double Individual::GetTarget()
{
	return target;
}
double Individual::GetIndTarget()
{
	return individualtarget;
}
void Individual::SetMaxPayment(double value)
{
	maxpayment = value;
	return;
}
double Individual::GetMaxPayment()
{
	return maxpayment;
}
void Individual::Mutator(int mp, int r, std::mt19937_64 & rng, double sigma)
{
	
	double temp=0;
    std::uniform_real_distribution<double> uniR(0.,1.);
    
    if(mp==0)
    {
        if(sigma<0){
            //   std::cout<<sigma<<" 1\n";
            genebelow[r]=uniR(rng)*maxpayment;
        }
        else   {
            
            std::normal_distribution<double> gaus(genebelow[r],sigma);
            //  std::cout<<sigma<<" 2\n";
            do
            {
                
                
                temp=gaus(rng);
            }while (((temp<0)||(temp>maxpayment)));
            genebelow[r]=temp;
        }
        
        
        if(r==0)
        {
            geneabove[r]=genebelow[r];
        }
    }
    if((mp==1)&&(r>0))//all have a threshold to mutate, but if it is the first round there is no point to change
    {
        
        std::normal_distribution<double> gaus(genethreshold1[r],sigma);
        do
        {temp=gaus(rng);
        }while (((temp<0)||(temp>1)));
        
        genethreshold1[r]=temp;
    }
    if(mp==2)//mutates the payment option
    {
        
        if(sigma<0) {geneabove[r]=uniR(rng)*maxpayment;}
        else  {
            
            std::normal_distribution<double> gaus(geneabove[r],sigma);
            do
            {temp=gaus(rng);
            }while (((temp<0)||(temp>maxpayment)));
            geneabove[r]=temp;
        }
        
        
        if(r==0)
        {
            genebelow[r]=geneabove[r];
        }
    }
	
	return;
	
}
double Individual::CalcPayment(int roundN, double const & totalPotsize)
{
	if(roundN==0)
	{return genebelow[roundN];}
	else if((genethreshold1[roundN]*target)<= totalPotsize)
	{return geneabove[roundN];}
	else 
	{return genebelow[roundN];}
}
