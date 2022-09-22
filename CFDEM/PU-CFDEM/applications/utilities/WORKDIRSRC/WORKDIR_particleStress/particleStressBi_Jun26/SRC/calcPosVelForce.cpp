#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <math.h>
#include "calcPosVelForce.h"

using namespace std;

void calcForce(particleCloud& particles, double **& fcoll, double **& ftan, double& k_n, double& k_t, double& gamma_n, double& gamma_t, double***& delta_t, double& mu_f, double& dt, 
				bool& tangential_history, double *&liq, bool& liquid_transfer, double **& liquidVol, double& surf_tension, double& fluid_visc, double **& fcap, double **& fvisc,
				bool**& first_touch )
{	
	const double pi = 3.1415926535897;
	
	double ddist = 0;
	
	double *n_vec = new double[3];
	double *t_vec = new double[3];
	
	double *dist   = new double[3];
	double delta_n;
	
	double *v_r = new double[3];
	double *v_n = new double[3];
	double *v_t = new double[3];
	
	double *w_r = new double[3];
	double *v_trel = new double[3];
			
	double *F_n = new double[3];
	double *F_t = new double[3];	
		
	double shrmag; double fs; double fmu;
	
	// Liquid bridge variable
	double tilde_h;
	double tilde_liquidVol; 
	double const_A;
	double const_B; 
	double const_C; 
	double Ca;
	double ksi;
	
	// Initialise
	for(int j=0 ; j < particles.numP(); j++)
	{
		for (int dir=0;dir<3;dir++)
		{	
			fcoll[j][dir] = 0;
			 ftan[j][dir] = 0;		
 			 fcap[j][dir] = 0;
 			fvisc[j][dir] = 0;				 			 

			if(!tangential_history)
			{
				for(int i = 0; i < particles.numP(); ++i) 
	 			{
	 	        	delta_t[i][j][dir] = 0;
	 	    	}			
			}	 
		}	
	}	
	
	for(int j=0 ; j < particles.numP(); j++)
	{		
		for(int i=j+1 ; i < particles.numP(); i++)
		{			 
			for (int dir=0;dir<3;dir++) 
				dist[dir] = particles.position(i)[dir]-particles.position(j)[dir];
			   
			ddist = sqrt( dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2] );
			
			for (int dir=0;dir<3;dir++) 
				n_vec[dir] = dist[dir] / ddist;
			
			delta_n = ddist - ( particles.radius(i) + particles.radius(j) );
			
			for (int dir=0;dir<3;dir++) 
				v_r[dir] = particles.velocity(i)[dir]-particles.velocity(j)[dir];

			for (int dir=0;dir<3;dir++) 
				v_n[dir] = ( v_r[0] * n_vec[0] + v_r[1] * n_vec[1] + v_r[2] * n_vec[2] ) * n_vec[dir]; 
			
			// If there is collision					
			if( delta_n <= 0 )
			{					
				cout << " During collision relative distance = " << delta_n << " between particle " << i << " & " << j << endl;
				for (int dir=0;dir<3;dir++)
					cout << " Velocity of particle " << i << " " << particles.velocity(i)[dir] << " & " << j << " " << particles.velocity(j)[dir] << endl;
				
				for (int dir=0;dir<3;dir++)
				{	 
					v_t[dir]  = v_r[dir] - v_n[dir];
					F_n[dir]  = - k_n * delta_n * n_vec[dir] - gamma_n * v_n[dir];
				}
				
				// Liquid tranfer: Cohesion & viscous forces here 
			    if (liquid_transfer)
			    {
					if(!first_touch[i][j])
					{
						first_touch[i][j] = true;
						first_touch[j][i] = true;
						calcLiqIn(liq, i, j, liquidVol);
					}
					
					// Capillary force model
					if ( liquidVol[i][j] > 0 )
					{	
						cout << " Volume of liquid bridge = " << liquidVol[i][j] << endl;
						tilde_h = delta_n / particles.radius(j);
						tilde_liquidVol = liquidVol[i][j] / ( particles.radius(j) * particles.radius(j) * particles.radius(j) );
						const_A = - 1.1 * pow( tilde_liquidVol, -0.53 );
						const_B = - 0.019 * log( tilde_liquidVol ) + 0.48;
						const_C = 0.0042 * log ( tilde_liquidVol ) + 0.078;
						for (int dir=0;dir<3;dir++)
							fcap[j][dir] = ( exp ( const_A * tilde_h + const_B ) + const_C ) * n_vec[dir] * pi * particles.radius(j) * surf_tension ;
					
						// Viscous force model
						Ca = fluid_visc * sqrt(   ( v_r[0] * n_vec[0] ) * ( v_r[0] * n_vec[0] )	
											    + ( v_r[1] * n_vec[1] ) * ( v_r[1] * n_vec[1] )	
											    + ( v_r[2] * n_vec[2] ) * ( v_r[2] * n_vec[2] ) ) / surf_tension;
					
						ksi = 1 - 1./sqrt( 1 + 2 * tilde_liquidVol / ( pi * tilde_h * tilde_h) );
					
						for (int dir=0;dir<3;dir++)
							fvisc[j][dir] = ( 3./2. * Ca / tilde_h * ksi * ksi ) * n_vec[dir] * pi * particles.radius(j) * surf_tension ;
						
						// Add to normal force
						for (int dir=0;dir<3;dir++)
							F_n[dir] -= fcap[j][dir]; // + fvisc[j][dir];	
					}		
				}	 
		
				for (int dir=0;dir<3;dir++)
					w_r[dir] = ( particles.radius(j) - 0.5 * ddist ) * particles.omega(j)[dir] / particles.radius(j);
				
				// Tangential relative velocity
				v_trel[0] = v_t[0] - ( dist[2]*w_r[1] - dist[1]*w_r[2] );	
				v_trel[1] = v_t[1] - ( dist[0]*w_r[2] - dist[2]*w_r[0] );
				v_trel[2] = v_t[2] - ( dist[1]*w_r[0] - dist[0]*w_r[1] );	

			    // shear history effects
			    if (tangential_history)
			    {
			        for (int dir=0;dir<3;dir++)
						delta_t[i][j][dir] += v_trel[dir] * dt;
			
			        // rotate shear displacements
					for (int dir=0;dir<3;dir++)					
			        	delta_t[i][j][dir] -= ( delta_t[i][j][0] * dist[0] + delta_t[i][j][1] * dist[1] + delta_t[i][j][2] * dist[2] ) 
							                * dist[dir] / particles.radius(j) / particles.radius(j) ;
			    }
		
			    shrmag = sqrt( delta_t[i][j][0] * delta_t[i][j][0] + delta_t[i][j][1] * delta_t[i][j][1] + delta_t[i][j][2] * delta_t[i][j][2]);

			    // tangential forces = shear + tangential velocity damping
				for (int dir=0;dir<3;dir++)
					F_t[dir] = - k_t * delta_t[i][j][dir]; 	

			    // rescale frictional displacements and forces if needed
			    fs = sqrt( F_t[0] * F_t[0] + F_t[1] * F_t[1] + F_t[2] * F_t[2] );
			    fmu = mu_f * sqrt( F_n[0] * F_n[0] + F_n[1] * F_n[1] + F_n[2] * F_n[2] );

			    // energy loss from sliding or damping
			    if (fs > fmu) {
			        if (shrmag != 0.0) 
					{
			            for (int dir=0;dir<3;dir++)
						{
							F_t[dir] *= fmu/fs;
			            	delta_t[i][j][dir] = - F_t[dir] / k_t;
			            }
			        }
			        else F_t[0] = F_t[1] = F_t[2] = 0.0;
			    }
			    else
			    {
			        for (int dir=0;dir<3;dir++)
						F_t[dir] = - gamma_t * v_t[dir];
			    }	
		
				for (int dir=0;dir<3;dir++)
				{	
					fcoll[j][dir] -= ( F_n[dir] + F_t[dir] );
					 ftan[j][dir] -= F_t[dir];
					 
					fcoll[i][dir] += ( F_n[dir] + F_t[dir] );
					 ftan[i][dir] += F_t[dir];						 
				}
			}
			else if( ( delta_n > 0 && liquidVol[i][j] > 0 ) ) //|| ( particles.liquidOn(i) > 0 ) || ( particles.liquidOn(j) > 0 ) )
			{
				cout << " During breaking up the bridge relative distance = " << delta_n << " between particle " << i << " & " << j << endl;
				for (int dir=0;dir<3;dir++)
					cout << " Velocity of particle " << i << " " << particles.velocity(i)[dir] << " & " << j << " " << particles.velocity(j)[dir] << endl;
				
				cout << " During breaking up volume of liquid bridge = " << liquidVol[i][j] << endl;
				
				double dist_rup = pow(liquidVol[i][j], 1./3.);
				if( delta_n < dist_rup )
				{
					// Capillary force model
					tilde_h = delta_n / particles.radius(j);
					tilde_liquidVol = liquidVol[i][j] / ( particles.radius(j) * particles.radius(j) * particles.radius(j) );
					const_A = - 1.1 * pow( tilde_liquidVol, -0.53 );
					const_B = - 0.019 * log( tilde_liquidVol ) + 0.48;
					const_C = 0.0042 * log ( tilde_liquidVol ) + 0.078;
					for (int dir=0;dir<3;dir++)
						fcap[j][dir] = ( exp ( const_A * tilde_h + const_B ) + const_C ) * n_vec[dir] * pi * particles.radius(j) * surf_tension ;
										
					// Viscous force model
					Ca = fluid_visc * sqrt(   ( v_r[0] * n_vec[0] ) * ( v_r[0] * n_vec[0] )	
										    + ( v_r[1] * n_vec[1] ) * ( v_r[1] * n_vec[1] )	
										    + ( v_r[2] * n_vec[2] ) * ( v_r[2] * n_vec[2] ) ) / surf_tension;
					
					ksi = 1 - 1./sqrt( 1 + 2 * tilde_liquidVol / ( pi * tilde_h * tilde_h) );
					
					for (int dir=0;dir<3;dir++)
					{	
						fvisc[j][dir] = ( 3./2. * Ca / tilde_h * ksi * ksi ) * n_vec[dir] * pi * particles.radius(j) * surf_tension ;
						fvisc[i][dir] = fvisc[j][dir] ; // For post-processing
					}
					// Add to normal force
				
					for (int dir=0;dir<3;dir++)
					{
						fcoll[j][dir] -= ( - fcap[j][dir] - fvisc[j][dir] );						
						fcoll[i][dir] += ( - fcap[j][dir] - fvisc[j][dir] );
					}
																
				}
				else
				{
					cout << " Test delta_n > dist_rup " << endl;
					calcLiqBack(liq, i, j, liquidVol);
					first_touch[i][j] = false;
					first_touch[j][i] = false;	
				}					
			}												
		}	
	}
}

// Calculate positions & velocities Sundar-110 (Multi-scale approach, pp.32)
void calcPosVel(particleCloud& particles, double**& pos, double**& vel, double**& omega, double**&fcoll, double**& ftan, double& dt)
{
	double mass_p;
	double I_p;
	const double pi = 3.1415926535897;	
	for(int i=0; i < particles.numP(); i++)
	{
		mass_p = 4./3. * pi * particles.radius(i) * particles.radius(i) * particles.radius(i) * particles.density(i);
		I_p = 2./5. * mass_p * particles.radius(i) * particles.radius(i);
		
		for(int j=0; j < 3; j++)
		{
			pos[i][j]   =   particles.position(i)[j] 
				          + particles.velocity(i)[j] * dt
					   	  +( 2./3. * particles.force(i)[j] / mass_p - 1./6. * fcoll[i][j] / mass_p ) * dt * dt;   
			
			vel[i][j]   =   particles.velocity(i)[j] 
				          + particles.force(i)[j] / mass_p * dt
					   	  +( 5./6. * particles.force(i)[j] / mass_p - 1./6. * fcoll[i][j] / mass_p ) * dt * dt;
	
			omega[i][j] =   particles.omega(i)[j] 
				          + particles.forceTan(i)[j] / I_p * dt
					   	  +( 5./6. * particles.forceTan(i)[j] / I_p - 1./6. * ftan[i][j] / I_p ) * dt * dt;
			
		}	
	}
}

// Calculate liquid volumes in liquid bridges
void calcLiqIn(double*& liq, int& i, int& j, double**& liquidVol)
{
	liquidVol[i][j] = liq[i] / 6. + liq[j] / 6.;
	liq[i] -= liq[i] / 6.;
	liq[j] -= liq[j] / 6.;	
}			

// Re-distribute liquid to particles
void calcLiqBack(double *&liq, int& i, int& j, double**& liquidVol)
{
	liq[i] += liquidVol[i][j]/2.;
	liq[j] += liquidVol[i][j]/2.;	
	liquidVol[i][j] = 0;	
}