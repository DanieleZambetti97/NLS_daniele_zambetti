/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __fluid_
#define __fluid_

//Random numbers generator
	#include "random.h"
	int seed[4];
	Random rnd;

//parameters, observables
	const int n_observable = 4;
	int iu, ih, im, ix;
	double walker[n_observable];

// averages
	double glob_av [n_observable],	glob_av2 [n_observable];
	double blk_av [n_observable],	blk_norm,	accepted,	attempted;
	double stima_u,	stima_h,	stima_m,	stima_x;
	double err_u,		err_h,	err_m,	err_x;

//configuration
	const int nspin = 50;
	double s[nspin];

// thermodynamical state
	double temp, j, h;

// simulation
int nstep, nblk, metro, restart, nRestart;

//functions
double Input(int, int, double, double, double, int, int);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(int);
void ConfFinal(void);
void Measure(void);
double Boltzmann(int, int);
int Pbc(int);
double Error(double,double,int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
