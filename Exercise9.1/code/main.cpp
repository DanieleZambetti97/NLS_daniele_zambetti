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
	
	double sum;
	ofstream Results("Best.dat");
	for(int l=0;l<nPopulation;l++){
		sum=0;
		for(int i=0;i<nElement;i=i+2) {
			int a = Selection();
			int b = Selection();
			BuildParents(a,b);
			Mutation( ParentsList[0] );
			Mutation( ParentsList[1] );
			Crossover(i);
		}

		for(int i=0;i<nElement;i++)	copy ( NewPopulation[i].begin(), NewPopulation[i].end(), Population[i].begin() );
		Check();
		Order();
		
		for(int i=0;i<nElement/2;i++)	{
			sum += TripList[i].length;
			//cout<<TripList[i].length<<"	"<<sum<<endl;
		}
		Results<<l<<"	"<<TripList[0].length<<"	"<<sum*2./double(nElement)<<endl;
/*		cout<<"The best of "<<l<<"-th population: "<<endl;
		Show(TripList[0],Population);
		cout<<endl;
*/	}

	Results.close();
	ofstream out("BestTrip.dat");
	for(int i=0;i<nCity;i++){
		out<<Population[TripList[0].name][i]<<"	"<<CityPosX[Population[TripList[0].name][i]]<<"	"<<CityPosY[Population[TripList[0].name][i]]<<endl;
	}
	out.close();
	out.open("BestLength.dat");
	out<<TripList[0].length<<endl;
	out.close();

return 0;
}



void Input (){	
	//leggo le posizioni delle città dal file CityPosition.dat
	ifstream CityPosition ("CityPosition.dat");

	for(int i=0; i<nCity; i++){
		CityPosition>>CityName[i];
		CityPosition>>CityPosX[i];
		CityPosition>>CityPosY[i];
		//cout<<"Nome della città: "<<CityName[i]<<" Posizione ("<<CityPosX[i]<<";"<<CityPosY[i]<<")"<<endl;
	}
	CityPosition.close();
	
	//creo la Population
	Population.resize(nElement);
	NewPopulation.resize(nElement);
	
	for(int i=0; i<nElement; i++){
		Population[i].resize(nCity);
		NewPopulation[i].resize(nCity);
		for(int j=0;j<nCity;j++){
			Population[i]=CityName;
			random_shuffle( Population[i].begin(), Population[i].end() );
		}
	}

	Parents.resize(2);
	for(int i=0;i<2;i++) Parents[i].resize(nCity);
	Scambia.resize(nCity);
}



void Order (){
	//creo il vettore di Trip
	for(int i=0; i<nElement; i++){
		TripList[i].name=i;
		TripList[i].length=Length(Population[i]);
	}
	sort (TripList.begin(), TripList.end());		//ordina TripList dal cammino più corto al più lungo
}

double Length (vector<int> v){
	double sum=0;
	int j;
	for(int i=0;i<nCity;i++){
		j=i+1;
		if(i==nCity-1){j=0;}
		sum += sqrt( pow(CityPosX[v[i]]-CityPosX[v[j]], 2.) + pow(CityPosY[v[i]]-CityPosY[v[j]], 2.) );
	}
	return sum;
}



void Check (){
	for(int j=0;j<nElement;j++) CheckFunction(TripList[j]);
}

int CheckFunction (Trip& x){
	int name=x.name;
	int i=0;
	vector<int> verifica = Population[name];
	sort (verifica.begin(), verifica.end());
	while(i<nCity){
		if(verifica[i]!=CityName[i]){
			cout<<"Il cammino "<<x.name<<" non rispetta i vincoli \n Trip:"<<endl;
			for(int j=0;j<nCity;j++) cout<<Population[name][j]<<" ";
			cout<<endl;
			for(int j=0;j<nCity;j++) cout<<verifica[j]<<" ";
			cout<<endl;
			break;
		}
		i++;
	}
	return 0;
}



int Selection(){
	return int( nElement*pow(rnd.Rannyu(),2.) );
}



int BuildParents(int a,int b){
	ParentsList[0]=TripList[a];
	ParentsList[1]=TripList[b];
	Parents[0]=Population[TripList[a].name];
	Parents[1]=Population[TripList[b].name];
	ParentsList[0].name=0;
	ParentsList[1].name=1;	
	return 0;
}



int Mutation(Trip& x ){
	
	if(rnd.Rannyu()<=0.1) Scambio ( Parents[x.name] );
	if(rnd.Rannyu()<=0.1) ScambioContiguo ( Parents[x.name] );
	if(rnd.Rannyu()<=0.1) Shift ( Parents[x.name] );
	if(rnd.Rannyu()<=0.1) ScambioParte ( Parents[x.name] );
	if(rnd.Rannyu()<=0.1) Inverti ( Parents[x.name] );
	return 0;}



int Crossover (int i){
	int cut = int ( rnd.Rannyu(0.,nCity) );
	Pairs.resize ( nCity-cut );
	Scambia=Parents[0];
	int t=i+2;
	double a=rnd.Rannyu();
	for(;i<t;i++){
		if(a<0.6){
			for (int j=0;j<nCity-cut;j++){
				Pairs[j].Vcity = Parents[0][cut+j];
				auto it = find ( Parents[1].begin(), Parents[1].end(), Parents[0][cut+j] );
				Pairs[j].icity = distance ( Parents[1].begin(), it );
			}
		
			sort( Pairs.begin(), Pairs.end() );
		}
		NewPopulation[i]=Parents[0];
		if(a<0.6){
			for(int j=0;j<nCity-cut;j++){
				NewPopulation[i][cut+j]=Pairs[j].Vcity;
			}
		}
		Parents[0]=Parents[1];
		Parents[1]=Scambia;
	}	
	return 0;
}



int Show (Trip& x, vector<vector<int>>& y){
		cout<<" Element: "<<x.name<<" Length: "<<x.length<<" Trip: ";
		for(int i=0; i<nCity; i++) cout<<y[x.name][i]<<" ";
		cout<<endl;
		return 0;
}

int Scambio (vector<int>& x){
	int a = int(rnd.Rannyu(0,nCity));
	int b = int(rnd.Rannyu(0,nCity));
	swap(x[a],x[b]);
	return 0;
}

int ScambioContiguo(vector<int>& x){
	int i = int(rnd.Rannyu(0,nCity-1));
	swap(x[i], x[i+1]);
	return 0;
}

int Shift (vector<int>& x){
	int shift = int(rnd.Rannyu(1,nCity));
	rotate(x.begin(),x.begin()+shift,x.end());
	return 0;
}

int ScambioParte (vector<int>& x){
	int inizio = int ( rnd.Rannyu(0,nCity/2) );
	int length = int( rnd.Rannyu(1,nCity/2-inizio) );
	int shift = int ( rnd.Rannyu(length,nCity/2-1) );
	vector<int> y(length);
	vector<int> z(length);
	for(int i=0;i<length;i++) z[i]=x[inizio+shift+i];
	for(int i=0;i<length;i++) y[i]=x[inizio+i];
	for(int i=0;i<length;i++) x[inizio+i]=z[i]; 
	for(int i=0;i<length;i++) x[inizio+shift+i]=y[i]; 
	return 0;
}

int Inverti (vector<int>& x){
	int inizio = int(rnd.Rannyu(0,nCity));
	int n = int(rnd.Rannyu(0,nCity-inizio));
	reverse(x.begin()+inizio,x.begin()+inizio+n);
	return 0;
}
