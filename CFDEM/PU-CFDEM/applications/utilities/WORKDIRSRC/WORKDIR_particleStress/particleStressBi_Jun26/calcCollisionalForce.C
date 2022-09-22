#include "calcCollisionalForce.H"

namespace Foam
{

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

	void calcForce(cfdemCloud& particles, 	int collisionModel,
						int j, 		// Particle	
						int i,  	// neighboring particles
				        	double **& fcoll, 
						double **& ftan,
						double& k_n, 
						double& k_t, 
						double& gamma_n, 
						double& gamma_t,
						double*& youngsModulus,		  
						double*& poissonsRatio,
						double**& coefficientRestitution,		
						double**& coefficientFriction,					 
						//double***& delta_t, 
						double*& delta_t, 
						double& mu_f, 
						const scalar& dt, 
						bool& tangential_history, 
						double *&liq, 
						bool& liquid_transfer, 
						double **& liquidVol, 
						double& surf_tension, 
						double& fluid_visc, 
						double **& fcap, 
						double **& fvisc,
						//bool**& first_touch,
						bool& first_touch,
						bool& cohesion,
						double& minimumDistanceVdW,
						double **& cohEnergyDens,
						double **& fcoh,
						const scalar& rhop,
						double& e_n,
						double& e_t,	   	    
					        symmTensor& sigma_coll_JI,
						const bool& bidisperse,
						const int & typeJ,
						const int & typeI,
						const bool & verbose 		)
	{	
		
		const double pi = constant::mathematical::pi; 
		
		double ddist(0);

		double *n_vec = new double[3];
		double *t_vec = new double[3];

		double *dist   = new double[3];
		double delta_n(0);

		double *v_r = new double[3];
		double *v_n = new double[3];
		double *v_t = new double[3];

		double *w_r = new double[3];
		double *v_trel = new double[3];

		double *F_n = new double[3];
		double *F_t = new double[3];	

		double shrmag(0); double fs(0); double fmu(0);

		// Liquid bridge variables
		double tilde_h(0);
		double tilde_liquidVol(0); 
		double const_A(0);
		double const_B(0); 
		double const_C(0); 
		double Ca(0);
		double ksi(0);

		// Cohesion variables
		double ri(0); double rj(0); 
		double smin(0); double s(0); double smax(0);

		double S_n(0); double S_t(0);
		double Beta(0); 
		double Y_star(0); double G_star(0); double R_star(0); double m_star(0);
        	
		e_n = coefficientRestitution[typeJ][typeI];		 	
		mu_f = coefficientFriction[typeJ][typeI];
				
        	m_star = 1. / ( 1./(4./3.*pi*rhop*pow(particles.radius(i),3)) + 1./(4./3.*pi*rhop*pow(particles.radius(j),3) ) );

		// Collisional stress
		sigma_coll_JI = symmTensor(0,0,0,0,0,0);
		
		if(collisionModel==0)
		{
			R_star = 1. / ( 1./particles.radius(i) + 1./particles.radius(j) );
			Y_star = 1. / ( (1-pow(poissonsRatio[typeI],2))/youngsModulus[typeI] + (1-pow(poissonsRatio[typeJ],2))/youngsModulus[typeJ] );
			G_star = 1. / ( 2*(2+poissonsRatio[typeI])*(1-poissonsRatio[typeI])/youngsModulus[typeI] + 2*(2+poissonsRatio[typeJ])*(1-poissonsRatio[typeJ])/youngsModulus[typeJ] );
			Beta = log(e_n)/sqrt(log(e_n)*log(e_n)+pi*pi);
		}
		
		

		// Initialise
		for (int dir=0;dir<3;dir++)
		{				 			 
			F_n[dir] = 0;
			F_t[dir] = 0;
			n_vec[dir] = 0;
			t_vec[dir] = 0;

			dist[dir] = 0;

			v_r[dir] = 0;
			v_n[dir] = 0;
			v_t[dir] = 0;

			w_r[dir] = 0;
			v_trel[dir] = 0;
			
			if(!tangential_history)
			{
				//delta_t[i][j][dir] = 0;
				delta_t[dir] = 0;
			}	 
		}
		
		for (int dir=0;dir<3;dir++) 
			dist[dir] = particles.position(i)[dir]-particles.position(j)[dir];	

		ddist = sqrt( dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2] );

		for (int dir=0;dir<3;dir++) 
			n_vec[dir] = dist[dir] / ( ddist );
								
		delta_n = ddist - ( particles.radius(i) + particles.radius(j) );

		for (int dir=0;dir<3;dir++) 
			v_r[dir] = particles.velocity(i)[dir]-particles.velocity(j)[dir];

		for (int dir=0;dir<3;dir++) 
			v_n[dir] = ( v_r[0] * n_vec[0] + v_r[1] * n_vec[1] + v_r[2] * n_vec[2] ) * n_vec[dir]; 	
			
		// If there is collision					
		if( delta_n <= 0 )
		{					
			if(verbose)
			{
				Info << " " << endl;
				Info << " During collision relative distance = " << delta_n << " between particle " << j << "(type "<<typeJ+1<<")" << " and " << i << "(type "<<typeI+1<<")" << endl;
				//for (int dir=0;dir<3;dir++)
					//Info << " Vel_" << dir << " of particle = " << i << " " << particles.position(i)[dir] << "\t and " << j << " " << particles.position(j)[dir] << endl;
			}
			
                	if(collisionModel==0)
                	{
				S_n = 2.*Y_star*sqrt(R_star*mag(delta_n));
				S_t = 8.*G_star*sqrt(R_star*mag(delta_n));
				k_n = 4./3.*Y_star*sqrt(R_star*mag(delta_n));
				k_t = 8*G_star*sqrt(R_star*mag(delta_n));
				gamma_n = -2.*sqrt(5./6.)*Beta*sqrt(S_n*m_star);				
				gamma_t = -2.*sqrt(5./6.)*Beta*sqrt(S_t*m_star);
				if(verbose)
				{
                       			Info << " k_n = " << k_n;
                        		Info << " gamma_n = " << gamma_n;
					Info << " k_t = " << k_t;
                        		Info << " gamma_t = " << gamma_t << endl;
				}	
			 }
                	else if (collisionModel==1)
                	{
				gamma_n *= m_star;
				gamma_t *= m_star;	 
                	}	

			for (int dir=0;dir<3;dir++)
			{	 
				v_t[dir]  = v_r[dir] - v_n[dir];
				F_n[dir]  = - k_n * delta_n * n_vec[dir] - gamma_n * v_n[dir];				
			}

			// Liquid tranfer: Cohesion & viscous forces here 
			if (liquid_transfer)
			{
				/*
				if(!first_touch[i][j])
				{
					first_touch[i][j] = true;
					first_touch[j][i] = true;
					calcLiqIn(liq, i, j, liquidVol);
				}
				*/
				if(!first_touch)
				{
					first_touch = true;
					calcLiqIn(liq, i, j, liquidVol);
				}				

				// Capillary force model
				if ( liquidVol[i][j] > 0 )
				{	
					if(verbose) Info << " Volume of liquid bridge = " << liquidVol[i][j] << endl;
					
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

			if(cohesion)
			{
				ri = particles.radius(i);
				rj = particles.radius(j);
				smin = minimumDistanceVdW;    

				for (int dir=0;dir<3;dir++)
				{
					fcoh[j][dir] = cohEnergyDens[typeJ][typeI]/3.0 
							* (2.0*ri*rj*(ri+rj+smin))/(smin*smin*(2*ri + 2*rj + smin)*(2*ri + 2*rj + smin)) 
							* ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1)   
							* ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1) * n_vec[dir] ;			
					F_n[dir] -= fcoh[j][dir];
				}
			}	 

			for (int dir=0;dir<3;dir++)
				w_r[dir] = ( particles.radius(j) - 0.5 * ddist ) * particles.omega(j)[dir] / particles.radius(j);

			// Tangential relative position
			v_trel[0] = v_t[0] - ( dist[2]*w_r[1] - dist[1]*w_r[2] );	
			v_trel[1] = v_t[1] - ( dist[0]*w_r[2] - dist[2]*w_r[0] );
			v_trel[2] = v_t[2] - ( dist[1]*w_r[0] - dist[0]*w_r[1] );	

			// shear history effects
			if (tangential_history)
			{
				for (int dir=0;dir<3;dir++)
						//delta_t[i][j][dir] += v_trel[dir] * dt;
						delta_t[dir] += v_trel[dir] * dt;
						
				// rotate shear displacements
					for (int dir=0;dir<3;dir++)					
					//delta_t[i][j][dir] -= ( delta_t[i][j][0] * dist[0] + delta_t[i][j][1] * dist[1] + delta_t[i][j][2] * dist[2] ) 
					//				* dist[dir] / particles.radius(j) / particles.radius(j) ;
					delta_t[dir] -= ( delta_t[0] * dist[0] + delta_t[1] * dist[1] + delta_t[2] * dist[2] ) 
									* dist[dir] / particles.radius(j) / particles.radius(j) ;

			}

			//shrmag = sqrt( delta_t[i][j][0] * delta_t[i][j][0] + delta_t[i][j][1] * delta_t[i][j][1] + delta_t[i][j][2] * delta_t[i][j][2]);
			shrmag = sqrt( delta_t[0] * delta_t[0] + delta_t[1] * delta_t[1] + delta_t[2] * delta_t[2]);

			// tangential forces = shear + tangential position damping
			for (int dir=0;dir<3;dir++)
				//F_t[dir] = - k_t * delta_t[i][j][dir]; 	
				F_t[dir] = - k_t * delta_t[dir]; 
					
			// rescale frictional displacements and forces if needed
			fs = sqrt( F_t[0] * F_t[0] + F_t[1] * F_t[1] + F_t[2] * F_t[2] );
			fmu = mu_f * sqrt( F_n[0] * F_n[0] + F_n[1] * F_n[1] + F_n[2] * F_n[2] );

			// energy loss from sliding or damping
			if (fs > fmu) 
			{
				if (shrmag != 0.0) 
				{
					for (int dir=0;dir<3;dir++)
					{
							F_t[dir] *= fmu/fs;
					       //delta_t[i][j][dir] = - F_t[dir] / k_t;
					             delta_t[dir] = - F_t[dir] / k_t;
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

				if(verbose)
				{
					Info << " F_n[" << dir << "]= " << setw(IOstream::defaultPrecision() + 6) << F_n[dir];
					Info << " F_t[" << dir << "]= " << setw(IOstream::defaultPrecision() + 6) << F_t[dir];
					Info << " v_j[" << dir << "]= " << setw(IOstream::defaultPrecision() + 6) << particles.velocity(j)[dir] << " " 
					     << " v_i[" << dir << "]= " << setw(IOstream::defaultPrecision() + 6) << particles.velocity(i)[dir];
					Info << " x_j[" << dir << "]= " << setw(IOstream::defaultPrecision() + 6) << particles.position(j)[dir] << " " 
					     << " x_i[" << dir << "]= " << setw(IOstream::defaultPrecision() + 6) << particles.position(i)[dir] << endl;
				}
				//fcoll[i][dir] += ( F_n[dir] + F_t[dir] );
				// ftan[i][dir] += F_t[dir];						 
			}
			
		}
		else if(delta_n > 0)
		{ 
			if(liquid_transfer) //|| ( particles.liquidOn(i) > 0 ) || ( particles.liquidOn(j) > 0 ) )
			{
				if(liquidVol[i][j] > 0 )
				{
					if(verbose) Info << " During breaking up the bridge relative distance= " << delta_n << " between particle " << i << "\t and " << j << endl;
					for (int dir=0;dir<3;dir++)
					{
						if(verbose) Info << " position of particle= " << i << " " << particles.position(i)[dir] << "\t and " << j << " " << particles.position(j)[dir] << endl;
					}
					if(verbose) Info << " During breaking up volume of liquid bridge = " << liquidVol[i][j] << endl;

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
							//fcoll[i][dir] += ( - fcap[j][dir] - fvisc[j][dir] );
						}

					}
					else
					{
						if(verbose) Info << " Test delta_n > dist_rup " << endl;
						calcLiqBack(liq, i, j, liquidVol);
						//first_touch[i][j] = false;
						//first_touch[j][i] = false;
						first_touch = true;	
					}
				}						
			}

			else if(cohesion)
			{
				ri = particles.radius(i);
				rj = particles.radius(j);
				smin = minimumDistanceVdW;  
				s = delta_n; // separating distance between the surfaces of the two interacting particles
				smax = (ri + rj) / 4.0; // To speed up the simulation, a maxmimum cutoff separation equal to d/4 is introduced. Beyond smax, the van der Waals cohesive force is not considered.	


				if (s < smax)
				{		
					if (s < smin)
					{	
        					for (int dir=0;dir<3;dir++)
						{
							fcoh[j][dir] =  cohEnergyDens[typeJ][typeI]/3.0 
									* (2.0*ri*rj*(ri+rj+smin))/(smin*smin*(2*ri + 2*rj + smin)*(2*ri + 2*rj + smin)) 
									* ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1)   
									* ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1) * n_vec[dir];
						      fcoll[j][dir] -= fcoh[j][dir];		
						}			
        				} 
        				else
					{ 	
        					for (int dir=0;dir<3;dir++)
						{
							fcoh[j][dir] =  cohEnergyDens[typeJ][typeI]/3.0 
									* (2.0*ri*rj*(ri+rj+s))/(s*s*(2*ri + 2*rj + s)*(2*ri + 2*rj + s)) 
									* ((s * (2*ri + 2*rj + s)) / ((ri + rj + s)*(ri + rj + s) - (ri-rj)*(ri-rj)) - 1)
									* ((s * (2*ri + 2*rj + s)) / ((ri + rj + s)*(ri + rj + s) - (ri-rj)*(ri-rj)) - 1) * n_vec[dir];
						      fcoll[j][dir] -= fcoh[j][dir];		
						}			

					}

				}		

			}
		
		}
		
		// Particle collisional stresses
		/*
		sigma_coll.xx() += ( particles.position(i)[0] - particles.position(j)[0] ) * ( F_n[0] + F_t[0] ); // fcoll[j][0]; 													
		sigma_coll.xy() += ( particles.position(i)[0] - particles.position(j)[0] ) * ( F_n[1] + F_t[1] ); // fcoll[j][1]; 													
		sigma_coll.xz() += ( particles.position(i)[0] - particles.position(j)[0] ) * ( F_n[2] + F_t[2] ); // fcoll[j][2]; 													
		sigma_coll.yy() += ( particles.position(i)[1] - particles.position(j)[1] ) * ( F_n[1] + F_t[1] ); // fcoll[j][1]; 													
		sigma_coll.yz() += ( particles.position(i)[1] - particles.position(j)[1] ) * ( F_n[2] + F_t[2] ); // fcoll[j][2];
		sigma_coll.zz() += ( particles.position(i)[2] - particles.position(j)[2] ) * ( F_n[2] + F_t[2] ); // fcoll[j][2];
		*/
				
		sigma_coll_JI.xx() = 0.5 * dist[0] * ( F_n[0] + F_t[0] ); // fcoll[j][0]; 																										
		sigma_coll_JI.yy() = 0.5 * dist[1] * ( F_n[1] + F_t[1] ); // fcoll[j][1]; 													
		sigma_coll_JI.zz() = 0.5 * dist[2] * ( F_n[2] + F_t[2] ); // fcoll[j][2];				
		sigma_coll_JI.xy() = 0.5 * dist[0] * ( F_n[1] + F_t[1] ); // fcoll[j][1]; 													
		sigma_coll_JI.xz() = 0.5 * dist[0] * ( F_n[2] + F_t[2] ); // fcoll[j][2]; 
		sigma_coll_JI.yz() = 0.5 * dist[1] * ( F_n[2] + F_t[2] ); // fcoll[j][2];
		
		if( delta_n <= 0 )
		{
			if(verbose)
			{
				Info << " " << endl;
				Info << " sigma_coll_JI.xx() = " << sigma_coll_JI.xx() << endl; 																										
				Info << " sigma_coll_JI.yy() = " << sigma_coll_JI.yy() << endl;													
				Info << " sigma_coll_JI.zz() = " << sigma_coll_JI.zz() << endl;				
				Info << " sigma_coll_JI.xy() = " << sigma_coll_JI.xy() << endl; 													
				Info << " sigma_coll_JI.xz() = " << sigma_coll_JI.xz() << endl;
				Info << " sigma_coll_JI.yz() = " << sigma_coll_JI.yz() << endl;	
			}	
		}
		
	
	}
	
}
