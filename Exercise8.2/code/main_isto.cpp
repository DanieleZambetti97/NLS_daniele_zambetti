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
#include <string>
#include <cmath>

#include "random.h"
#include "metropolis.h"
#include "datablocking.h"

using namespace std;

double Psi(double, double, double);
double Psi_Second( double, double, double);
double Potential(double);

int main (int argc, char *argv[]){
	if(argc !=3){
		cout<<" error: lanciare il comando ./mani.exe mu sigma normalization";
	}

	int M=100000;
	int N=1000;
	int L=M/N;
	double AV[L];

	double Start=0.0;
	double ene,sum, mu, sigma;
	mu = atof(argv[1]);
	sigma = atof(argv[2]);
	Metropolis Met(Start,mu,sigma,1.);
	Met.Equilibration(0.50);

	ofstream pos("Posizioni.dat");
	for (int l=0;l<L;l++){
		sum=0;
		for (int i=0;i<N; i++){
			Met.MakePass();
			//calcolo l'energia in questo nuovo punto
			ene = -0.5*Psi_Second(Met.GetPosition(),mu,sigma)/Psi(Met.GetPosition(),mu,sigma)+Potential(Met.GetPosition());
			sum += ene;
			//salvo le posizioni
			pos << Met.GetPosition() <<endl;
		}
		AV[l]=sum/double(N);
	}
	pos.close();

	Data_Blocking (AV, L);

//stampo la stima di Epot ottenuta dal datablocking su Energie.dat
	ifstream risul("Risultati.dat");
	ofstream energie("Energie.dat",ios::app);
	double x;
	while(! risul.eof()){
		risul>>x;
	}
	energie<<mu<<"	"<<sigma<<"	"<<x<<endl;
	energie.close();
	risul.close();
	return 0;
}


double Psi(double x,double mu, double sigma){
	return (exp( -(x-mu)*(x-mu)/(2.*sigma*sigma) ) + exp( -(x+mu)*(x+mu)/(2.*sigma*sigma) ) );
}
double Psi_Second(double x, double mu, double sigma){
	return (exp( -(x-mu)*(x-mu)/(2.*sigma*sigma) )*( pow((x-mu)/(sigma*sigma),2.) - 1./pow(sigma,2.) ) + exp( -(x+mu)*(x+mu)/(2.*sigma*sigma) ) * ( pow((x+mu)/(sigma*sigma),2.) - 1./pow(sigma,2.) ) );
}
double Potential(double x){
	return pow(x,4.)-5./2.*pow(x,2.);
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
