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
#include <cstring>
#include <cmath>

#include "random.h"
#include "datablocking.h"

using namespace std;

Random rnd=initialization();
double T=1;
double K=100;
double mu=0.1;
double sigma=0.25;


double GBM_Direct ( double, double, double, double);
double GBM_Step ( double, double, double, double, double);
double Max ( double, double);

int main (int argc, char *argv[]){

	int M=10000;
	int N=100;
	int L=M/N;
	double AV[L];
	double sum;

	string outputfile;
	double S_t, FinalPrice;

	for(int Graph=0;Graph<4;Graph++){		//cicla sui 4 grafici richiesti (put-discretized, call-discretized, put-direct, call-direct)

		for (int j=0; j<L; j++){				//cicla sul n di blocchi
			sum = 0;
			for(int i=0; i<N; i++){				//cicla sulle stime del prezzo che in ogni blocco vengono fatte

				if(Graph<2){
					S_t = K;
					for(double h=0; h<1; h = h+0.01)  S_t = GBM_Step(h, h+0.01, mu, sigma, S_t);
				}

				FinalPrice = exp(-mu*T) * Max(Graph, S_t);

				sum+=FinalPrice;
			}

			AV[j]=sum/double(N);
		}
	outputfile = "Risultati"+to_string(Graph+1)+".dat";
	const char * c = outputfile.c_str();
	DataBlocking( AV, L, c);
	}
	return 0;
}

double GBM_Direct (double T, double mu, double sigma, double K){
	return K* exp( (mu-1./2.*sigma*sigma)*T + sigma*rnd.Gauss(0.,T) );
}

double GBM_Step (double T_0, double T_1, double mu, double sigma, double S_Prima){
	return S_Prima* exp( (mu-1./2.*sigma*sigma)*(T_1-T_0) + sigma*rnd.Gauss(0.,1.)*sqrt(T_1-T_0) );
}
double Max (double i, double S_t){
	if(i==0) return max(0. , K-S_t);
	if(i==1) return max(0. , S_t-K);
	if(i==2) return max(0. , K-GBM_Direct(T,mu,sigma,K));
	if(i==3) return max(0. , GBM_Direct(T,mu,sigma,K)-K);
	else return 0;
};
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
