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
#include "datablocking.h"

using namespace std;

int main (int argc, char *argv[]){

   Random rnd=initialization();
	int M=100000;
	int N=1000;
	int L=M/N;
	double PI[L];
	int Nhit=0;
	int Ntot=0;
	double y_1, y_2, theta;
	double l=0.8;
	double t=1.;

	for (int j=0; j<L; j++){

		Nhit=Ntot=0;

		for(int i=0; i<N; i++){

			theta=rnd.Direction2D();
			y_1=rnd.Rannyu()*10.;
			y_2=y_1+l*sin(theta);

			if (int(y_1)!=int(y_2) || y_1<0 || y_2<0) {Ntot++;Nhit++;}
			else {Ntot++;}
		}

		PI[j]=2.*l*Ntot/(Nhit*t);
	}

	char outputfile[]="Risultati.dat";
	DataBlocking(PI, L, outputfile);	

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
