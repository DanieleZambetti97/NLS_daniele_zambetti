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

#include "random.h"

using namespace std;

Random rnd=initialization();

double NewPos(double i, double* pos);

int main (int argc, char *argv[]){

	int nWalk=10000;
	int nStep=100;
	double AVDist[nStep]={0};
	double AVDist2[nStep]={0};
	double ErrorAVDist[nStep];
	double pos[3];
	double dist;

	for(int i=0;i<2;i++){
		for(int l=0;l<nWalk;l++){

			fill(&pos[0],&pos[3],0.);
			for(int k=0;k<nStep;k++){
				NewPos(i, pos);
				dist=pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2];
				AVDist[k] += dist;
				AVDist2[k] += dist*dist;
			}
		}

		ofstream risul ("Risultati"+to_string(i)+".dat");
		for(int i=0;i<nStep;i++){
			AVDist[i] = AVDist[i]/double(nWalk);   //trasformo le somme in medie
			AVDist2[i] = AVDist2[i]/double(nWalk);
			ErrorAVDist[i] = sqrt( (AVDist2[i] - AVDist[i]*AVDist[i])/double(nWalk) );
			risul << i+1 << "	" << sqrt(AVDist[i]) << "	" << ErrorAVDist[i] / (2.*sqrt(AVDist[i])) << endl;
		}

		risul.close();
	}


	return 0;
}

double NewPos(double i, double* pos){

	double x,y,z;

	if(i==0){
		x = int(rnd.Rannyu()*6.);
		if(x==0){pos[0]+= 1.;}
		if(x==1){pos[0]+=-1.;}
		if(x==2){pos[1]+= 1.;}
		if(x==3){pos[1]+=-1.;}
		if(x==4){pos[2]+= 1.;}
		if(x==5){pos[2]+=-1.;}
	}
	if(i==1){
		x=y=z=1.;
		while(x*x+y*y+z*z>1){
			x = rnd.Rannyu(-1,1);
			y = rnd.Rannyu(-1,1);
			z = rnd.Rannyu(-1,1);
		}
		x = x/sqrt(x*x+y*y+z*z);
		y = y/sqrt(x*x+y*y+z*z);
		z = z/sqrt(x*x+y*y+z*z);
		pos[0] += x;
		pos[1] += y;
		pos[2] += z;
	}
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
