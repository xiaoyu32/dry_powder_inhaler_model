#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <math.h>

//#include "particleCloud.h"
#include "calcPosVelForce.h"

using namespace std;
using std::scientific;

int main(int argc, char* argv[])
{	
	// Read simulations parameters
    string char_endTime;
    string char_outputFreq;
	double endTime;
	double outputFreq;
	ifstream sim_props("sim_props.dat");
    if (sim_props.is_open())
    {
		sim_props >> char_endTime >> endTime; 
		sim_props >> char_outputFreq >> outputFreq;	
	  	sim_props.close();
    }
    else cout << "Unable to open file sim_props.dat"; 
	
	// Initialise variables	
	int MaxNumberOfParticles = 10; 
	#include "init_variables.h"
	
	// Particle properties & initial positions & velocities
	string Number_of_particles;
	string header;
    int numberOfParticles;
	ifstream parts_init ("parts_init.dat");
    if (parts_init.is_open())
    {
		parts_init >> Number_of_particles;
		parts_init >> numberOfParticles;
		cout << Number_of_particles << " " << numberOfParticles << endl;
		parts_init >> header;
		for(int line = 0; line < numberOfParticles; ++line)
		{
		    parts_init >> pos[line][0] >> pos[line][1] >> pos[line][2] 
				       >> vel[line][0] >> vel[line][1] >> vel[line][2] 
				       >> omega[line][0] >> omega[line][1] >> omega[line][2] 
					   >> rad[line] 
					   >> liq[line] 
					   >> density[line];			
		}	
	  	parts_init.close();
    }

    else cout << "Unable to open file parts_init,dat"; 

	// Create particle cloud
	particleCloud particles;
	particles.setNumberOfParticles (numberOfParticles);
	
	// Effective mass
	const double pi = 3.1415926535897;	
	for(int i=0 ; i < particles.numP(); i++)
	{	
		mass[i] = density[i] * 4./3. * pi * rad[i] * rad[i] * rad[i];
	}	
	
	// Read collision parameters
	double k_n;
	double k_t;
	double gamma_n;
	double gamma_t;
	double mu_f;
	double e_n;
	double e_t;
	bool tangential_history;
	string str_tangential_history;
	bool liquid_transfer;
	string str_liquid_transfer;
	double N_coll; 
	double surf_tension;
	double fluid_visc;

	string char_coll;
    ifstream coll_props("coll_props.dat");
    if (coll_props.is_open())
    {
		coll_props >> char_coll >> k_n;
			cout << char_coll << " = " << k_n << endl;
		coll_props >> char_coll >> e_n;
			cout << char_coll << " = " << e_n << endl;
		coll_props >> char_coll >> e_t;
			cout << char_coll << " = " << e_t << endl;	
		coll_props >> char_coll >> mu_f;
			cout << char_coll << " = " << mu_f << endl;			
		coll_props >> char_coll >> N_coll;
			cout << char_coll << " = " << N_coll << endl;			
		coll_props >> char_coll >> str_tangential_history;	
		if( str_tangential_history == "false" )
		{
			tangential_history = false;
			cout << "Tangential shear history is off " << endl;
		}else
		{
			tangential_history = true;
			cout << "Tangential shear history is on " << endl;
		}	
		coll_props >> char_coll >> str_liquid_transfer;	
		if( str_liquid_transfer == "false" )
		{
			liquid_transfer = false;
			cout << "Liquid transfer is off " << endl;
		}else
		{
			liquid_transfer = true;
			cout << "Liquid transfer is on " << endl;
		}		
		coll_props >> char_coll >> surf_tension;
			cout << char_coll << " = " << surf_tension << endl;			
		coll_props >> char_coll >> fluid_visc;
			cout << char_coll << " = " << fluid_visc << endl; 
      	coll_props.close();
    }
	
    else cout << "Unable to open file coll_props.dat"; 	
		
	// Allocate variables
	particles.allocateVariables();
	
	// Time-step
	double m_eff = 0;
	double dt = 10.e8;
	double dt_sim = 0;
	double t_coll = 0; 
	double min_t_coll = 0;
	
	// Calculate time step due to collision time 
	for(int j=0 ; j < particles.numP(); j++)
	{	
		for(int i=j+1 ; i < particles.numP(); i++)
		{
			// --> effective mass into m_eff[1,2]
			m_eff = ( mass[i]*mass[j] ) / ( mass[i] + mass[j] );
			k_t = 2./7. * k_n;
			gamma_n = - 2 * log(e_n) * sqrt( m_eff * k_n / ( pi * pi + log(e_n) * log(e_n) ) );
			gamma_t = - 2 * log(e_t) * sqrt( m_eff * k_t / ( pi * pi + log(e_t) * log(e_t) ) );
			t_coll = sqrt ( ( pi * pi + log(e_t) * log(e_t) ) / ( k_n / m_eff ) ) ;
				if ( min_t_coll < t_coll ) min_t_coll = t_coll;	 		
			dt_sim = t_coll / N_coll;
				if ( dt_sim < dt ) dt = dt_sim;	 
		}	
	}	
	
	cout << "k_t = " << k_t << endl;
	cout << "gamma_n = " << gamma_n << endl;
	cout << "gamma_t = " << gamma_t << endl;
	
	cout << "Collision time = " << min_t_coll << " Simulation time-step = " << dt << endl;

	// Create output files
	ofstream particleFile_1("particle_1.dat");		
	particleFile_1 << "Time[s] " << '\t' << "Position_x " << '\t' << "Position_y " << '\t' << "Position_z " 
								 << '\t' << "Velocity_x " << '\t' << "Velocity_y " << '\t' << "Velocity_z "
							     << '\t' << "Omega_x    " << '\t' << "Omega_y    " << '\t' << "Omega_z    " 
								 << '\t' << "F_coll_x   " << '\t' << "F_coll_y   " << '\t' << "F_coll_z   " 
								 << '\t' << "F_tan_x    " << '\t' << "F_tan_y    " << '\t' << "F_tan_z    " 
							     << '\t' << "Radius     " 
								 << '\t' << "Density    " 
								 << '\t' << "LiquidOn   " 
								 << '\t' << "F_cap_x    " << '\t' << "F_cap_y    " << '\t' << "F_cap_z    "
								 << '\t' << "F_visc_x   " << '\t' << "F_visc_y   " << '\t' << "F_visc_z   "
								 << endl;	

	ofstream particleFile_2("particle_2.dat");		
	particleFile_2 << "Time[s] " << '\t' << "Position_x " << '\t' << "Position_y " << '\t' << "Position_z " 
								 << '\t' << "Velocity_x " << '\t' << "Velocity_y " << '\t' << "Velocity_z "
							     << '\t' << "Omega_x    " << '\t' << "Omega_y    " << '\t' << "Omega_z    " 
								 << '\t' << "F_coll_x   " << '\t' << "F_coll_y   " << '\t' << "F_coll_z   " 
								 << '\t' << "F_tan_x    " << '\t' << "F_tan_y    " << '\t' << "F_tan_z    " 
							     << '\t' << "Radius     " 
								 << '\t' << "Density    " 
								 << '\t' << "LiquidOn   " 
								 << '\t' << "F_cap_x    " << '\t' << "F_cap_y    " << '\t' << "F_cap_z    "
								 << '\t' << "F_visc_x   " << '\t' << "F_visc_y   " << '\t' << "F_visc_z   "
								 << endl;		
	
	// Number of timestep
	int Ndt = endTime/dt;
	double simTime;
	
	// Time integration 
	for(int Idt=1; Idt < Ndt; Idt++ )
	{
		// Simulation time
		simTime = Idt * dt;
		cout << " Time[s] = " << simTime << endl;
		
		particles.setPos(pos);
		particles.setVel(vel, omega);		
		particles.setRadius(rad);
		particles.setDensity(density);
								 	
		// Calculate collisional force
		calcForce(			  particles, 
								  fcoll, 
								   ftan, 
								    k_n, 
									k_t, 
								gamma_n, 
								gamma_t, 
								delta_t, 
								   mu_f, 
								     dt, 
					 tangential_history, 
							  		liq, 
						liquid_transfer, 
						      liquidVol, 
						   surf_tension, 
						     fluid_visc, 
							       fcap, 
								  fvisc,
						    first_touch 	);
					   
		particles.setForce(fcoll, ftan, fcap, fvisc);
		
		// Set liquid on particle
		particles.setLiquid(liq);		
		
		// Calculate positions & velocities
		calcPosVel(particles, pos, vel, omega, fcoll, ftan, dt);

		for(int i=0 ; i < particles.numP(); i++)
		{
			if( i == 0)
			{	
				particleFile_1 << std::scientific 
				 << simTime 	 << " "			 
				 << particles.position(i)[0] << " " << particles.position(i)[1] << " " << particles.position(i)[2] << '\t'
				 << particles.velocity(i)[0] << " " << particles.velocity(i)[1] << " " << particles.velocity(i)[2] << '\t'
				 << particles.omega(i)[0]    << " " << particles.omega(i)[1]    << " " << particles.omega(i)[2]    << '\t'
				 << particles.force(i)[0]    << " " << particles.force(i)[1]    << " " << particles.force(i)[2]    << '\t'
				 << particles.forceTan(i)[0] << " " << particles.forceTan(i)[1] << " " << particles.forceTan(i)[2] << '\t'
				 << particles.radius(i)	  	 << '\t'
				 << particles.density(i)	 << '\t' 					 
				 << particles.liquidOn(i)    << '\t'
				 << particles.forceCap(i)[0] << " " << particles.forceCap(i)[1] << " " << particles.forceCap(i)[2] << '\t'
				 << particles.forceVis(i)[0] << " " << particles.forceVis(i)[1] << " " << particles.forceVis(i)[2] 					 
				 << " "	 				     << endl ;					
			}
			else
			{
				particleFile_2 << std::scientific 
				 << simTime    <<	" "			 
				 << particles.position(i)[0] << " " << particles.position(i)[1] << " " << particles.position(i)[2] << '\t'
				 << particles.velocity(i)[0] << " " << particles.velocity(i)[1] << " " << particles.velocity(i)[2] << '\t'
				 << particles.omega(i)[0]    << " " << particles.omega(i)[1]    << " " << particles.omega(i)[2]    << '\t'
				 << particles.force(i)[0]    << " " << particles.force(i)[1]    << " " << particles.force(i)[2]    << '\t'
				 << particles.forceTan(i)[0] << " " << particles.forceTan(i)[1] << " " << particles.forceTan(i)[2] << '\t'					 
				 << particles.radius(i)	  	 << '\t'
				 << particles.density(i)	 << '\t' 					 
				 << particles.liquidOn(i)    << '\t'
				 << particles.forceCap(i)[0] << " " << particles.forceCap(i)[1] << " " << particles.forceCap(i)[2] << '\t'
				 << particles.forceVis(i)[0] << " " << particles.forceVis(i)[1] << " " << particles.forceVis(i)[2] 					 
				 << " "	 				     << endl ;				
			}
		}		
	}
	particleFile_1.close();
	particleFile_2.close();
	
				   					      
	cout << " :-) " << endl ;
}