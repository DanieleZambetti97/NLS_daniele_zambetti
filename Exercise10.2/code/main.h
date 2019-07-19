
#ifndef __main__
#define __main__

#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

//variabili da inizializare
bool circle=false;
int NCity=30;
double beta=1;
int NStep=1000;
int NTemp=2000;
int size,myrank;
Random rnd;

vector<double> CityPosX(NCity);
vector<double> CityPosY(NCity);
vector<int>	CityName(NCity);
vector<int> Path;
vector<int> proposed_path;


//funzioni utilizzate nel programma
void Input();
int CheckFunction(vector<int>&);
int Mutation(vector<int>&);

double Length (vector<int>&);
int Scambio (vector<int>&);
int ScambioContiguo(vector<int>&);
int Shift (vector<int>&);
int ScambioParte (vector<int>&);
int Inverti (vector<int>&);

#endif 
