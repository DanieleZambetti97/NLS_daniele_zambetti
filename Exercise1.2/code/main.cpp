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
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd=initialization();
	
	int M=10000;
	int N;
	int n[4]={1,2,10,100};
	double SumExp=0;
	double SumLorentz=0;
	double SumGauss=0;

	ofstream exp("Risultati_Exp.dat");
	ofstream lor("Risultati_Lor.dat");
	ofstream gauss("Risultati_Gauss.dat");

	for(int k=0;k<4;k++){
		N=n[k];
		for(int j=0;j<M;j++){	
			SumExp=0;
			SumLorentz=0;
			SumGauss=0;
			for(int i=0;i<N;i++){
				SumExp+=rnd.Exp(1.);
				SumLorentz+=rnd.Lorentz(1.,0.);
				SumGauss+=rnd.Gauss(1.,0.5);
				}
			exp<<SumExp/double(N)<<endl;
			lor<<SumLorentz/double(N)<<endl;
			gauss<<SumGauss/double(N)<<endl;
		}
	}

	exp.close();
	lor.close();
	gauss.close();
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
