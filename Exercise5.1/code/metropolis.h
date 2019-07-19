/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Metropolis__
#define __Metropolis__

#include "random.h"

class Metropolis {
	private:
		double Pos[3];
		double delta;
		double Acceptance;
		int Ntot; 
		int Nacc;
		int f;
		Random rnd=initialization();

	protected:
	public:
	// constructors
	Metropolis (double*,double, int);
	// destructor
	~Metropolis();
	// methods
	bool MakePass();
	double Equilibration(double);		
	double* GetPosition();
	double GetAcceptance();
};

double PD_1s (double,double,double);
double PD_2p (double,double,double);
//bool Metropolis(double);


#endif // __Metropolis__


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
