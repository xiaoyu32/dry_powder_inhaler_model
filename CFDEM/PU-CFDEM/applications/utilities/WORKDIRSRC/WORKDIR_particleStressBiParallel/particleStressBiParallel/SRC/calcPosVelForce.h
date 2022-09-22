#include "particleCloud.h"
void calcPosVel(particleCloud&, double**&, double**&, double**&, double**&, double**&, double&);
void calcForce(    particleCloud&, //particles, 
                   double **&, // fcoll, 
				   double **&, // ftan, 
				   double&   , // k_n, 
				   double&   , // k_t, 
				   double&   , // gamma_n, 
				   double&   , // gamma_t, 
				   double***&, // delta_t, 
				   double&   , // mu_f, 
				   double&   , // dt, 
				   bool&     , // tangential_history, 
				   double  *&, // liq, 
				   bool&     , // liquid_transfer, 
				   double **&, // liquidVol, 
				   double&   , // surf_tension, 
				   double&   , // fluid_visc, 
				   double **&, // fcap, 
				   double **&, // fvisc
				   bool   **&  // liquid bridge boolean 			
			   );
void calcLiqIn  (double*&, int&, int&, double **&);
void calcLiqBack(double*&, int&, int&, double **&);