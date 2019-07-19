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
#include "random.h"
#include "datablocking.h"

using namespace std;

double Data_Blocking ( double * Average, int Block_Number){
	ofstream risul("Risultati.dat");
	double finalAV[Block_Number];
	double finalAV2[Block_Number];
	double err[Block_Number];
	double ave2,sum;

	for (int j=0; j<Block_Number; j++){
		sum=ave2=0;
		for (int i=0; i<j+1; i++){
			sum+=Average[i];
			ave2+=Average[i]*Average[i];
		}
		finalAV[j]=sum/double(j+1);
		finalAV2[j]=ave2/double(j+1);
		if(j!=0){err[j]=error(finalAV[j],finalAV2[j],j);}
		else{err[j]=0;}
		risul<<double(j+1)<<"	"<<err[j]<<"	"<<finalAV[j]<<endl;
	}
	risul.close();
	return 0.;
}

double error(double SumAve, double SumAve2, int n){
	return sqrt( ( SumAve2-(SumAve*SumAve) )/double(n) );
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
