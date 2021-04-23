/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
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

#include "metropolis.h"
#include "random.h"

using namespace std;

Metropolis :: Metropolis (double* Position, double Delta,int i){
	for(int i=0;i<3;i++){
		Pos[i]=Position[i];
	}
	delta=Delta;
	Acceptance=0.;
	Ntot=0;
	Nacc=0;
	f=i;
}

Metropolis ::	~Metropolis() {}

bool Metropolis ::  MakePass(){
	double New_Pos[3];
	double A=0;
	for(int i=0;i<3;i++){
//		New_Pos[i]=rnd.Rannyu(Pos[i]-delta/2.,Pos[i]+delta/2.);
		New_Pos[i]=rnd.Gauss(Pos[i],delta);
	}
	if(f==0) {A=min(1., PD_1s(New_Pos[0],New_Pos[1],New_Pos[2])/PD_1s(Pos[0],Pos[1],Pos[2]) );}
	if(f==1) {A=min(1., PD_2p(New_Pos[0],New_Pos[1],New_Pos[2])/PD_2p(Pos[0],Pos[1],Pos[2]) );}
	double r=rnd.Rannyu();
	if(r<=A) {
		for(int i=0;i<3;i++){Pos[i]=New_Pos[i];}
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

double* Metropolis :: GetPosition(){
	return Pos;
}

double Metropolis :: GetAcceptance(){
	return Acceptance;
}

double Metropolis :: Equilibration(double Desire_Acceptance){
	cout<<"Equilibration: \nDesire acceptance = "<<Desire_Acceptance<<endl;
	cout<<" Acc 	  Delta  Numero_correzioni"<<endl;
	for(int i=0;i<10000;i++){
		if(Acceptance==Desire_Acceptance){
			cout<<Acceptance<<"	"<<delta<<"		"<<i<<endl;
			break;
		}
		Ntot=Nacc=0;
		for(int j=0;j<100;j++){
			MakePass();
		}
//		cout<<Acceptance<<"	"<<delta<<"		"<<i<<endl;
		if(Acceptance<Desire_Acceptance){delta=delta-1./(i+1.);}
		if(Acceptance>Desire_Acceptance){delta=delta+1./(i+1.);}
		if(i==9999){
			cout<<Acceptance<<"	"<<delta<<"		"<<i<<endl;
			cout<<"Dopo "<<i+1<<" tentativi non ho raggiunto il valore di accettazione desiderato \n Acceptance ="<<Acceptance<<endl;}
	}
	return Acceptance;
}

double PD_1s (double x,double y,double z){
	double r=sqrt(x*x+y*y+z*z);
	return M_PI*exp(-2*r);
}
double PD_2p (double x,double y,double z){
	double r=sqrt(x*x+y*y+z*z);
	double CosTheta=z/r;
	return 1./(32.*M_PI)*exp(-r)*pow(r*CosTheta,2);
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
