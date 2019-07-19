#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

//variabili da inizializare
bool circle=false;
int nCity=30;
int nElement=900;		//DEVE ESSERE PARI!!!!!!
int nPopulation=100;
Random rnd=initialization();


//struc utlizzate nel programma
struct Trip{
	int name;		//linea della matrice popolazione
	double length;	//lunghezza del percorso
	friend bool operator<(const Trip& l, const Trip& r){
        return l.length<r.length;
	}
};

struct Pair{
	int Vcity;
	int icity;
	friend bool operator<(const Pair& l, const Pair& r){
        return l.icity<r.icity;
	}
};


//funzioni utilizzate nel programma
vector<double> CityPosX(nCity);
vector<double> CityPosY(nCity);
vector<int>	CityName(nCity);
vector<vector<int>> Population;
vector<vector<int>> NewPopulation;
vector<vector<int>> Parents;
vector<Trip> TripList(nCity*nCity);
vector<Trip> ParentsList(2);
vector<Pair> Pairs;
vector<int> Scambia;

void Input();
void Order();
void Check();
int BuildParents(int,int);
int CopyParents(int);
int Selection();
int Mutation(Trip&);
int Crossover(int);

double Length (vector<int>);
int CheckFunction(Trip&);
int Show(Trip&,vector<vector<int>>&);
bool CrossoverOrder(int,int);
int Scambio (vector<int>&);
int ScambioContiguo(vector<int>&);
int Shift (vector<int>&);
int ScambioParte (vector<int>&);
int Inverti (vector<int>&);


