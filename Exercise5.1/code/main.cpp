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

#include "datablocking.h"
#include "random.h"
#include "metropolis.h"

using namespace std;

int main (int argc, char *argv[]){
	ofstream position("Posizioni.dat");
	int M=1000000;
	int N=10000;
	int L=M/N;
	double AV[L];
	double Start[3]={1.,1.,1.};
	double sum;

	double* Pos;

	for(int i=0;i<2;i++){
		Metropolis Met(Start,1.,i);
		Met.Equilibration(0.50);

		for (int l=0;l<L;l++){
			sum=0;
			for (int n=0;n<N;n++){
				Met.MakePass();
				Pos=Met.GetPosition();
				sum+= sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1]+Pos[2]*Pos[2]);
				if(((n+1)+(l+1)*N)%100 == 0){position<<Pos[0]<<"	"<<Pos[1]<<"	"<<Pos[2]<<"	"<<endl;}				//stampo 1 punto ogni 100
			}
			AV[l]=sum/double(N);
//			cout<<Met.GetAcceptance()<<endl;
		}
		string str="Risultati"+to_string(i)+".dat";
		const char * outputfile = str.c_str();
		DataBlocking(AV, L, outputfile);
	}
	position.close();

	return 0;
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
