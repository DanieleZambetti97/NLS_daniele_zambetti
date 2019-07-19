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
 
int main (int argc, char *argv[]){

	Random rnd = initialization();

	int M=10000;
	int N=100;
	int L=M/N;
	double AV[L];
	double sum;

	for (int j=0; j<L; j++){
		sum=0;
		for(int i=0; i<N; i++){
			sum += M_PI/2.*cos( M_PI/2.* rnd.Rannyu() );
		}
		AV[j] = sum/double(N);
	}

	char outputfile[] = "RisultatiUnifrome.dat";
	DataBlocking(AV, L, outputfile);


	double x,f_x;

	for (int j=0; j<L; j++){
		sum=0;
		for(int i=0; i<N; i++){
			x = 1.-sqrt (1.-rnd.Rannyu());							//genero un x con distribuzione p(x)=2-2x
			sum += M_PI/2.*cos( M_PI/2.*x)/(2*(1-x));
		}
		AV[j] = sum/double(N);
	}

	strcpy(outputfile, "RisultatiOttimiza.dat");
	DataBlocking(AV, L, outputfile);

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
