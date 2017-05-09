#include "spectdg.h"

// user input parameters:
#define EPS_i 1.0e-5 //relative error in spatial integration
#define JMAX_i 8 //maximum power of 2 in spatial integration grid
#define JMIN_i 8 //minimum power of 2 in spatial integration grid

#define pi_i 3.14159
#define omega_i 7.2921e-05 //omega=2*pi/(23*3600+56*60+4) earth's angular velocity
#define rho_i 1009.765 //average density
#define f1_i -0.4102 //initial phi
#define t1_i pi_i/6 //initial theta
#define t2_i 1.1262 //final theta
#define h_i 4.0e3 //ocean layer total depth
#define r1_i 6.370e6 //eath's radius at the bottom of ocean layer
#define M_i 3 //spectral terms in theta direction in large-scale
#define N_i 3 //spectral terms in phi direction in large-scale
#define M1_i 3 //spectral terms in theta direction in small-scale
#define N1_i 3 //spectral terms in phi direction in small-scale
#define dt0_i 2.0 //time step in terms of days
#define day_i 5000 //total time in days
#define ah_i 1.0e2 //horizontal viscosity coefficient
#define av_i -2.25e-8 //vertical viscosity coefficient

SpectDG::Input parameters()
// creates the input data structure based on user defined parameters and returns the structure
// returns a data structure of type SpectDG::Input
{
	SpectDG::Input input = {EPS_i,JMAX_i,JMIN_i,pi_i,omega_i,rho_i,f1_i,t1_i,t2_i,h_i,r1_i,M_i,N_i,M1_i,N1_i,dt0_i,day_i,ah_i,av_i};
	return input;
}