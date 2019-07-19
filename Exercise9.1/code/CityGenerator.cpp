#include <iostream>
#include <fstream>
#include <cmath>

#include "random.h"
#include "main.h"

using namespace std;

int main(){
ofstream Position ("CityPosition.dat");
Random rnd=initialization();

	if(circle==true){
		double x;
		for(int i=0; i<nCity; i++){
			x=rnd.Rannyu(0,2.*M_PI);
			Position<<i<<" ";
			Position<<sin(x) <<" ";
			Position<<cos(x) <<"\n";
		}
	}else{
		for(int i=0; i<nCity; i++){
			Position<<i<<" ";
			Position<<rnd.Rannyu()<<" ";
			Position<<rnd.Rannyu()<<"\n";
		}
	}
Position.close();

return 0;
}
