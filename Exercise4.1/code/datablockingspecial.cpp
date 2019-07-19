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
#include <string>
#include <cstdlib>

#include "datablockingspecial.h"

using namespace std;

double DataBlocking ( char * File_Input, int Dati_Number, char * File_Output){
	double Dati[Dati_Number];
	double AV[Dati_Number];
	double AV2[Dati_Number];
	double Err[Dati_Number];
	double ave2,sum;

	int k=0;double x;
	ifstream dati(File_Input);
	while(dati>>x) {
		k++;
	} 
	dati.close();
	dati.open(File_Input);
	for(int i=0;i<k;i++){
		dati>>x;
		if(i>=k-Dati_Number){
			Dati[i-k+Dati_Number]=x;
		}
	}
	dati.close();
	ofstream risul(File_Output);
	for (int j=0; j<Dati_Number; j++){
		sum=ave2=0;
		for (int i=0; i<j+1; i++){
			sum+=Dati[i];
			ave2+=Dati[i]*Dati[i];
		}
		AV[j]=sum/double(j+1);
		AV2[j]=ave2/double(j+1);
		Err[j]=error(AV[j],AV2[j],j);

		if(j==0) { Err[j]=0; }

		risul<<double(j+1)<<" "<<AV[j]<<" "<<Err[j]<<endl;
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
