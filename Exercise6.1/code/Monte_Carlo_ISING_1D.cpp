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
#include <ostream>
#include <cstdlib>
#include <cmath>

#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(int argc, char** argv) {
	if(argc!=8) cout << "Passa dalla riga di comando ./Monte_Carlo_ISING_1D.exe  Restart(0 or 1) Metropolis(0 or 1) Temperature  j h BlockNumber StepNumber" << endl;

	Input(atoi(argv[1]), atoi(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), atoi(argv[6]), atoi(argv[7])); //Inizialization

	for(int iblk=1; iblk <= nblk; ++iblk) { //Simulation
		Reset(iblk);   //Reset block averages
		for(int istep=1; istep <= nstep; ++istep) {
			Move(metro);
			Measure();
			Accumulate(); //Update block averages
		}
		Averages(iblk);   //Print results for current block
	}
	ConfFinal(); //Write final configuration

return 0;
}


double Input(int Restart, int  Metro, double Temp, double J, double H, int  nBlock, int  nStep) {

	restart = Restart;
	temp = Temp;
	metro = Metro;
	j = J;
	h = H;
	nblk = nBlock;
	nstep = nStep;

	ifstream ReadInput;
	if (restart==0) {
		nRestart=0;
	} else {
		ReadInput.open("Restart.dat");
		ReadInput >> nRestart;
		ReadInput.close();
	}

	cout << "Classic 1D Ising model             " << endl;
	cout << "Monte Carlo simulation             " << endl << endl;
	cout << "Nearest neighbour interaction      " << endl << endl;
	cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
	cout << "The program uses k_B=1 and mu_B=1 units " << endl;
	cout << "Temperature = " << temp << endl;
	cout << "Number of spins = " << nspin << endl;
	cout << "Exchange interaction = " << j << endl;
	cout << "External field = " << h <<endl;

	if(metro==1) cout << "The program perform Metropolis moves" << endl;
	else cout << "The program perform Gibbs moves" << endl;
	
	cout << "Number of blocks = " << nblk << endl;
	cout << "Number of steps in one block = " << nstep << endl << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();    

//Prepare arrays for measurements
	iu = 0; //Energy
	ih = 1; //Heat capacity
	im = 2; //Magnetization
	ix = 3; //Magnetic susceptibility
 

//initial configuration
	if(restart==0){
		for (int i=0; i<nspin; ++i) {
			if(rnd.Rannyu() >= 0.5) s[i] = 1;
			else s[i] = -1;
		}
	} else {
		ReadInput.open("config.final");
		for (int i=0; i<nspin; ++i) {ReadInput>>s[i];}
		ReadInput.close();
	}
	for(int i=0;i<nspin;i++){
//		cout<<s[i]<<endl;
	}
  
//Evaluate energy etc. of the initial configuration
	Measure();

//Print initial values for the potential energy and virial
	cout << "Initial energy = " << walker[iu]/(double)nspin << endl << endl;
	cout << "This is the " << nRestart+1 << " run" << endl << endl;

return 0;
}


void Move(int metro) {
	int o;
	double p, energy_old, energy_new;

	for(int i=0; i<nspin; ++i) {
		o = (int)(rnd.Rannyu()*nspin);	//Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)

 		if(metro==1){ //Metropolis
			energy_old = Boltzmann(s[o], o);
			s[o] = s[o]*-1.;
			energy_new = Boltzmann(s[o], o);
			p=exp(-(energy_new-energy_old) / temp);
			if( rnd.Rannyu()<=p ) { //accetiamo la mossa
				accepted++;
				attempted++;
			}else{
				s[o]=s[o]*-1.;	
				attempted++;
			}
		}else{ //Gibbs sampling
			p = exp(-Boltzmann(1., o )/temp) / (exp(-Boltzmann(1.,o)/temp)+ exp(-Boltzmann(-1.,o)/temp) );//probabilità che s[o] = +1
			if( rnd.Rannyu() <= p ) {
				s[o]=+1;
			}else{
				s[o]=-1;
			}
		}
	}
}


double Boltzmann(int sm, int ip) {
	double ene = -j * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
	return ene;
}

void Measure() {
	double u = 0.0, m=0.0;

	for (int i=0; i<nspin; ++i) {
		u += -j * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
		m += s[i];
	}
	walker[iu] = u;
	walker[ih] = u*u;
	walker[im] = m;
	walker[ix] = m*m;
}


void Reset(int iblk) { //Reset block averages
	if(iblk == 1) {
		for(int i=0; i<n_observable; ++i) {
			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}

	for(int i=0; i<n_observable; ++i) {
		blk_av[i] = 0;
	}
	blk_norm = 0;
	attempted = 0;
	accepted = 0;
}


void Accumulate(void) { //Update block averages

	for(int i=0; i<n_observable; ++i) {
		blk_av[i] = blk_av[i] + walker[i];
   }
	blk_norm = blk_norm + 1.0;
}


void Averages(int iblk){ //Print results for current block

	ofstream Ene, Heat, Mag, Chi;
    
	cout << "Block number " << iblk << endl;
	cout << "Acceptance rate " << accepted/attempted << endl << endl;

	//Energy
	Ene.open("./risultati/DataBlocking_Ene.dat",ios::app);
	stima_u = blk_av[iu]/blk_norm/double(nspin);
	glob_av[iu]  += stima_u;
	glob_av2[iu] += stima_u*stima_u;
	err_u=Error(glob_av[iu],glob_av2[iu],iblk);
	Ene << nRestart*nblk+iblk << "	" << glob_av[iu]/double(iblk) << "	" << err_u << endl;
	Ene.close();

	//Heat Capacity
	Heat.open("./risultati/DataBlocking_Heat.dat",ios::app);
	stima_h = 1./(temp*temp) * ( blk_av[ih]/blk_norm - (blk_av[iu]*blk_av[iu]) / (blk_norm*blk_norm) ) / double(nspin);
	glob_av[ih]  += stima_h;
	glob_av2[ih] += stima_h*stima_h;
	err_h=Error(glob_av[ih],glob_av2[ih],iblk);
	Heat << nRestart*nblk+iblk << "	" << glob_av[ih]/double(iblk) << "	" << err_h << endl;
	Heat.close();

	// Magnetic Susceptibility
	Chi.open("./risultati/DataBlocking_Susce.dat",ios::app);
	stima_x = 1./temp * blk_av[ix]/blk_norm / double(nspin); 
	glob_av[ix]  += stima_x;
	glob_av2[ix] += stima_x*stima_x;
	err_x = Error( glob_av[ix], glob_av2[ix], iblk);
	Chi << nRestart*nblk+iblk << "	" << glob_av[ix]/double(iblk) << "	" << err_x << endl;
	Chi.close();

	// Magnetization
	Mag.open("./risultati/DataBlocking_Magn.dat",ios::app);
	stima_m = blk_av[im]/blk_norm / double(nspin); 
	glob_av[im]  += stima_m;
	glob_av2[im] += stima_m*stima_m;
	err_m = Error( glob_av[im], glob_av2[im], iblk);
	Mag << nRestart*nblk+iblk << "	" << glob_av[im]/double(iblk) << "	" << err_m << endl;
	Mag.close();

	cout << "----------------------------" << endl << endl;
}


void ConfFinal(void) {
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i) {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();
	
	nRestart++;
  WriteConf.open("Restart.dat");
	WriteConf << nRestart;
  rnd.SaveSeed();
}

int Pbc(int i) { //Algorithm for periodic boundary conditions
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk) {
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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
