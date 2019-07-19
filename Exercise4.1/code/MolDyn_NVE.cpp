/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <cstdlib>     // srand, rand: to generate random number
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <string>
#include <iomanip>

#include "MolDyn_NVE.h"
#include "datablockingspecial.h"

using namespace std;

int main(int argc, char* argv[]){
	
	if(argc !=2){
		cout<<"Start a new simulation 		= 0 \nRestart from last configuration = 1 \n : ";
		cin>>restart;
		if (restart==true) system("rm risultati/*");
	}else {
		restart=atoi(argv[1]);
	}


	Input(); //Inizialization
	int nconf = 1;
	for(int istep=1; istep <= nstep; ++istep){
		Move(); //Move particles with Verlet algorithm
		Measure(); //Properties measurement
		if(istep%100 == 0){
			cout << "Number of time-steps: " << istep << endl;
			//ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
			nconf += 1;
		}
	}
	ConfFinal();         //Write final configuration to restart

	DataBlocking ("./risultati/output_epot.dat", nstep, "./risultati/ave_epot.dat");
	DataBlocking ("./risultati/output_etot.dat", nstep, "./risultati/ave_etot.dat");
	DataBlocking ("./risultati/output_ekin.dat", nstep, "./risultati/ave_ekin.dat");
	DataBlocking ("./risultati/output_temp.dat", nstep, "./risultati/ave_temp.dat");
	DataBlocking ("./risultati/output_pres.dat", nstep, "./risultati/ave_pres.dat");

	return 0;
}

void Input(void){ //Prepare all stuff for the simulation
	ifstream ReadInput,ReadConf;

	cout << "Classic Lennard-Jones fluid" << endl;
	cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "The program uses Lennard-Jones units " << endl;

	seed = 1;    //Set seed for random numbers
	srand(seed); //Initialize random number generator

	ReadInput.open(input_file); //Read input
	
	ReadInput >> temp;

	ReadInput >> npart;
	cout << "Number of particles = " << npart << endl;

	ReadInput >> rho;
	cout << "Density of particles = " << rho << endl;
	vol = (double)npart/rho;
	cout << "Volume of the simulation box = " << vol << endl;
	box = pow(vol,1.0/3.0);
	cout << "Edge of the simulation box = " << box << endl;

	ReadInput >> rcut;
	ReadInput >> delta;
	ReadInput >> nstep;
	ReadInput >> blocklength;	

	cout << "The program integrates Newton equations with the Verlet method " << endl;
	cout << "Time step = " << delta << endl;
	cout << "Number of steps = " << nstep << endl << endl;
	ReadInput.close();

	if(restart==false){ 
	//Read initial configuration
	  cout << "Read initial configuration from file config.0 " << endl << endl;
	  ReadConf.open("config.0");
	  for (int i=0; i<npart; ++i){
	    ReadConf >> x[i] >> y[i] >> z[i];
	    x[i] = x[i] * box;
	    y[i] = y[i] * box;
	    z[i] = z[i] * box;
	  }
	  ReadConf.close();
	
	  //Prepare initial velocities
	  cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
	  double sumv[3] = {0.0, 0.0, 0.0};
	  for (int i=0; i<npart; ++i){
	    vx[i] = rand()/double(RAND_MAX) - 0.5;
	    vy[i] = rand()/double(RAND_MAX) - 0.5;
	    vz[i] = rand()/double(RAND_MAX) - 0.5;

	    sumv[0] += vx[i];
	    sumv[1] += vy[i];
	    sumv[2] += vz[i];
	  }
	  for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
	  double sumv2 = 0.0, fs;
	  for (int i=0; i<npart; ++i){
	    vx[i] = vx[i] - sumv[0];
	    vy[i] = vy[i] - sumv[1];
	    vz[i] = vz[i] - sumv[2];

	    sumv2 += vx[i]*	vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	  }
	  sumv2 /= (double)npart;
	
	  fs =sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
	  for (int i=0; i<npart; ++i){
	    vx[i] *= fs;
	    vy[i] *= fs;
	    vz[i] *= fs;
	
	    xold[i] = Pbc(x[i] - vx[i] * delta);
	    yold[i] = Pbc(y[i] - vy[i] * delta);
	    zold[i] = Pbc(z[i] - vz[i] * delta);
	  }
	}
	if(restart==true){
	//Read final configuration and restart from there
		ifstream ReadFinalConf,ReadOldFinalConf;
		cout << "Restart form final configuration from file config.final " << endl << endl;
		ReadFinalConf.open("config.final");
		ReadOldFinalConf.open("config.old_final");
		double startx[npart], starty[npart], startz[npart];
		for (int i=0; i<npart; ++i){
			ReadFinalConf >> x[i] >> y[i] >> z[i];
			ReadOldFinalConf >> xold[i] >> yold[i] >> zold[i];
			x[i] = x[i] * box;
			y[i] = y[i] * box;
			z[i] = z[i] * box;
			xold[i] = xold[i] * box;
			yold[i] = yold[i] * box;
			zold[i] = zold[i] * box;
			startx[i] = x[i];
			starty[i] = y[i];
			startz[i] = z[i];
		}
		ReadFinalConf.close();
		ReadOldFinalConf.close();
		
	//Prepare initial velocities
		double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];
		double fs, sumv2=0;

		for(int i=0; i<npart; ++i){ //Force acting on particle i
			fx[i] = Force(i,0);
			fy[i] = Force(i,1);
			fz[i] = Force(i,2);
		}

		for(int i=0; i<npart; ++i){ //calcolo v(t+dt/2)

			xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
			ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
			znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

			vx[i] = Pbc(xnew - x[i])/ delta;
			vy[i] = Pbc(ynew - y[i])/ delta;
			vz[i] = Pbc(znew - z[i])/ delta;

			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];	

			x[i] = xnew;
			y[i] = ynew;
			z[i] = znew;
		}
		sumv2 /= (double)npart;
		fs = sqrt(3.*temp/sumv2);   // fs = velocity scale factor 
		
		for (int i=0; i<npart; ++i){
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;

			xold[i] = Pbc(startx[i] - vx[i] * delta);
			yold[i] = Pbc(starty[i] - vy[i] * delta);
			zold[i] = Pbc(startz[i] - vz[i] * delta);

			x[i] = startx[i];
			y[i] = starty[i];
			z[i] = startz[i];
   	}
	}
	
	//per fare datablocking su tutte le configurazioni (aggiunte in append ai fil gia esistenti) devo passargli un nconf che tiene conto di quante volte sono ripartito
	return;
}


void Move(void){ //Move particles with Verlet algorithm
	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

	for(int i=0; i<npart; ++i){ //Force acting on particle i
		fx[i] = Force(i,0);
		fy[i] = Force(i,1);
		fz[i] = Force(i,2);
	}

	for(int i=0; i<npart; ++i){ //Verlet integration scheme

		xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
		ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
		znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

		vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
		vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
		vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

		xold[i] = x[i];
		yold[i] = y[i];
		zold[i] = z[i];

		x[i] = xnew;
		y[i] = ynew;
		z[i] = znew;
	}
	return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
	double f=0.0;
	double dvec[3], dr;

	for (int i=0; i<npart; ++i){
		if(i != ip){
			dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
			dvec[1] = Pbc( y[ip] - y[i] );
			dvec[2] = Pbc( z[ip] - z[i] );

			dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
			dr = sqrt(dr);

			if(dr < rcut){
				f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
			}
		}
	}

	return f;
}

void Measure(){ //Properties measurement
	double v, t, p, vij, pij;
	double dx, dy, dz, dr;
	ofstream Epot, Ekin, Etot, Temp, Pres;

	Epot.open("./risultati/output_epot.dat",ios::app);
	Ekin.open("./risultati/output_ekin.dat",ios::app);
	Temp.open("./risultati/output_temp.dat",ios::app);
	Etot.open("./risultati/output_etot.dat",ios::app);
	Pres.open("./risultati/output_pres.dat",ios::app);

	v = 0.0; //reset observables
	t = 0.0;
	p = 0.0;

	//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i){
		for (int j=i+1; j<npart; ++j){
	
			dx = Pbc( x[i] - x[j] );
			dy = Pbc( y[i] - y[j] );
			dz = Pbc( z[i] - z[j] );

			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);
			pij= 1./pow(dr,12)-0.5/pow(dr,6);
			p += pij;

			if(dr < rcut){
				vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
				//Potential energy
				v += vij;
			}
		}          
	}

	for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]); //Kinetic energy
   
		stima_epot = v/(double)npart; //Potential energy
		stima_ekin = t/(double)npart; //Kinetic energy
		stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
		stima_etot = (t+v)/(double)npart; //Total enery
		stima_pres = rho * stima_temp + 48.* p / (3.*vol); //Pressure

		Epot << stima_epot << endl;
		Ekin << stima_ekin << endl;
		Temp << stima_temp << endl;
		Etot << stima_etot << endl;
		Pres << stima_pres << endl;

		Epot.close();
		Ekin.close();
		Temp.close();
		Etot.close();
		Pres.close();

	return;
}


void ConfFinal(void){ //Write final configuration
	ofstream WriteFinalConf;
	ofstream WriteOldFinalConf;
	cout << "Print final configuration to file config.final and scond to last configuration to file config.old_final " << endl << endl;
	WriteFinalConf.open("config.final");
	WriteOldFinalConf.open("config.old_final");

	for (int i=0; i<npart; ++i){
		WriteFinalConf << x[i]/box << "	" << y[i]/box << "	" << z[i]/box /*<< setw(wd) << vx[i] << setw(wd) <<  vy[i] << setw(wd) << vz[i]*/ << endl;
		WriteOldFinalConf << xold[i]/box << "	" <<  yold[i]/box << "	" << zold[i]/box << endl;
	}
	WriteFinalConf.close();
	WriteOldFinalConf.close();

	return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
	ofstream WriteXYZ;

	WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
	WriteXYZ << npart << endl;
	WriteXYZ << "This is only a comment!" << endl;
	for (int i=0; i<npart; ++i){
		WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
	}
	WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
	return r - box * rint(r/box);
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
