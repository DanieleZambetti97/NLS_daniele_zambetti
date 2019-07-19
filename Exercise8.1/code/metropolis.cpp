/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Si	mulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "metropolis.h"

using namespace std;

Metropolis :: Metropolis (double Position, double mu, double sigma, double delta){
	Start_Pos=Position;
	Pos=Position;
	Mu=mu;
	Sigma=sigma;
	Delta=delta;
	Acceptance=0.;
	Ntot=0;
	Nacc=0;
}

Metropolis ::	~Metropolis() {}

bool Metropolis ::  MakePass(){
	double New_Pos=rnd.Rannyu(Pos-Delta/2.,Pos+Delta/2.);
	double A=min(1., Psi2(New_Pos)/Psi2(Pos) );
	double r=rnd.Rannyu();
	if(r<=A) {
		Pos=New_Pos;
		Ntot++;
		Nacc++;
		Acceptance=double(Nacc)/double(Ntot);
		return true;
	}else{
		Ntot++;
		Acceptance=double(Nacc)/double(Ntot);
		return false;
	}
}

double Metropolis :: GetPosition(){
	return Pos;
}

double Metropolis :: GetAcceptance(){
	return Acceptance;
}

double Metropolis :: Equilibration(double Desire_Acceptance){
//	cout<<"Equilibration: \nDesire acceptance = "<<Desire_Acceptance<<endl;
//	cout<<" Acc 	  Delta  Numero_correzioni"<<endl;
	int i;
	for(i=0;i<1000;i++){
		if(Acceptance==Desire_Acceptance){break;}
		Ntot=Nacc=0;
		for(int j=0;j<100;j++){
			MakePass();
		}
//		cout<<Acceptance<<"	"<<delta<<"		"<<i<<endl;
		if(Acceptance<Desire_Acceptance){Delta=Delta-1./(i+1.);}
		if(Acceptance>Desire_Acceptance){Delta=Delta+1./(i+1.);}
		if(i==999){cout<<"Dopo "<<i<<" tentativi non ho raggiunto il valore di accettazione desiderato \n Acceptance = "<<Acceptance<<endl;}
	}
//	cout<<"N di tentativi: "<<i<<"\naccettazione raggiunta: "<<Acceptance<<endl;
	return Acceptance;
}

double Metropolis :: Psi2(double x){
	double a=exp(-pow((x-Mu)/Sigma,2.)*0.5);
	double b=exp(-pow((x+Mu)/Sigma,2.)*0.5);
	return (a*a+b*b+2.*a*b);
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
