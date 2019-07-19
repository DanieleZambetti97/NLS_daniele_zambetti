#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "random.h"
#include "main.h"

using namespace std;

int main () {
	Input();
	
	ofstream Results("Best.dat");
	for(int l=0;l<NTemp;l++){
		beta *= 1.003;
		int i=0;
		for(;i<NStep;i++){
			proposed_path=Path;
			Mutation(proposed_path);
			if(exp( -beta * (Length(proposed_path)-Length(Path)) )>=rnd.Rannyu()){
				Path=proposed_path;
				CheckFunction(Path);
			}
		}
		Results << beta << "	" << Length(Path) << endl;
	}
	Results.close();

	//stampo il percorso ottenuto
	Results.close();
	ofstream out("BestTrip.dat");
	for(int i=0;i<NCity;i++){
		out << Path[i] << "	" << CityPosX[Path[i]] << "	" << CityPosY[Path[i]] << endl;
	}
	out.close();
	out.close();
	out.open("BestLength.dat");
	out<<Length(Path)<<endl;
	out.close();

	return 0;
}


void Input (){	
	//leggo le posizioni delle città dal file CityPosition.dat
	ifstream CityPosition ("CityPosition.dat");

	for(int i=0; i<NCity; i++){
		CityPosition>>CityName[i];
		CityPosition>>CityPosX[i];
		CityPosition>>CityPosY[i];
		//cout<<"Nome della città: "<<CityName[i]<<" Posizione ("<<CityPosX[i]<<";"<<CityPosY[i]<<")"<<endl;
	}
	CityPosition.close();
	
	//creo il cammino iniziale
	Path=CityName;
	random_shuffle( Path.begin(), Path.end() );
}


double Length (vector<int>& v){
	double sum=0;
	int j;
	for(int i=0;i<NCity;i++){
		j=i+1;
		if(i==NCity-1){j=0;}
		sum += sqrt( pow(CityPosX[v[i]]-CityPosX[v[j]], 2.) + pow(CityPosY[v[i]]-CityPosY[v[j]], 2.) );
	}
	return sum;
}



int CheckFunction (vector<int>& x){
	int i=0;
	vector<int> verifica = x;
	sort (verifica.begin(), verifica.end());
	while(i<NCity){
		if(verifica[i]!=CityName[i]){
			cout<<"Il cammino non rispetta i vincoli \n Trip:"<<endl;
			for(int j=0;j<NCity;j++) cout<<x[i]<<" ";
			cout<<endl;
			for(int j=0;j<NCity;j++) cout<<verifica[j]<<" ";
			cout<<endl;
			break;
		}
		i++;
	}
	return 0;
}


int Mutation( vector<int>& x ){
	double r= rnd.Rannyu(0,5);
	if(r<=1) Scambio ( x );
	if(r>1 and r<=2) ScambioContiguo ( x );
	if(r>2 and r<=3) Shift ( x );
	if(r>3 and r<=4) ScambioParte ( x );
	if(r>4) Inverti ( x );
	return 0;
}



int Scambio (vector<int>& x){
	int a = int(rnd.Rannyu(0,NCity));
	int b = int(rnd.Rannyu(0,NCity));
	swap(x[a],x[b]);
	return 0;
}

int ScambioContiguo(vector<int>& x){
	int i = int(rnd.Rannyu(0,NCity-1));
	swap(x[i], x[i+1]);
	return 0;
}

int Shift (vector<int>& x){
	int shift = int(rnd.Rannyu(1,NCity));
	rotate(x.begin(),x.begin()+shift,x.end());
	return 0;
}

int ScambioParte (vector<int>& x){
	int inizio = int ( rnd.Rannyu(0,NCity/2) );
	int length = int( rnd.Rannyu(1,NCity/2-inizio) );
	int shift = int ( rnd.Rannyu(length,NCity/2-1) );
	vector<int> y(length);
	vector<int> z(length);
	for(int i=0;i<length;i++) z[i]=x[inizio+shift+i];
	for(int i=0;i<length;i++) y[i]=x[inizio+i];
	for(int i=0;i<length;i++) x[inizio+i]=z[i]; 
	for(int i=0;i<length;i++) x[inizio+shift+i]=y[i]; 
	return 0;
}

int Inverti (vector<int>& x){
	int inizio = int(rnd.Rannyu(0,NCity));
	int n = int(rnd.Rannyu(0,NCity-inizio));
	reverse(x.begin()+inizio,x.begin()+inizio+n);
	return 0;
}
