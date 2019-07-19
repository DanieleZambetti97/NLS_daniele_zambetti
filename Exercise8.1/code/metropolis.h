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

class Metropolis {
	private:
		double Start_Pos;
		double Pos;
		double Mu;
		double Sigma;
		double Delta;
		double Acceptance;
		int Ntot; 
		int Nacc;
		Random rnd=initialization();

	protected:
	public:
	// constructors
	Metropolis (double,double,double,double);
	// destructor
	~Metropolis();
	// methods
	bool MakePass();
	double Equilibration(double);		
	double GetPosition();
	double GetAcceptance();

	double Psi2(double); 
};
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
