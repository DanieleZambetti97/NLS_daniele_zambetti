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
#include <cstring>
#include <algorithm>

#include "random.h"
#include "datablocking.h"

using namespace std;

int main () {
	Random rnd = initialization();
	int M = 10000;
	int N = 100;
	int L = M/N;
	double AV[L];
	double sum;

	for (int j=0; j<L; j++){
		sum = 0;

		for(int i=0; i<N; i++){
			sum += rnd.Rannyu();
		}

		AV[j]=sum/double(N);
	}

	char OutputFile[] = "Risultati1.dat";
	DataBlocking( AV, L, OutputFile );
 
	for (int j=0; j<L; j++){
		sum=0;
		for(int i=0; i<N; i++){
			double a=rnd.Rannyu();
			sum+=(a-0.5)*(a-0.5);
		}
		AV[j]=sum/double(N);
	}

	strcpy( OutputFile, "Risultati2.dat" );
	DataBlocking( AV, L, OutputFile );

	int nBin=100;
	int nNumbers=10000;
	int Bin [nBin];
	ofstream risultati("RisultatiChiQuadro.dat");

	for(int k=0;k<100; k++){
		fill( Bin, Bin+nBin, 0);
		sum=0;
		for (int j=0; j<nNumbers; j++){
			int b=int( rnd.Rannyu()*nBin );
			Bin[ b ]++;
		}
		for(int i=0;i<nBin;i++){
			double a=double(nNumbers)/double(nBin);
			sum +=( (Bin[i]-a) * (Bin[i]-a) / a );
		}
		risultati<<k<<"	"<<sum<<endl;
	}
	risultati.close();
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
