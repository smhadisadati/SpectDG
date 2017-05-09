#include "spectdg.h"

void SpectDG::initialize(SpectDG::Input input)
// sets the input parameters for a SpectDG object and initializes the solution parameters
// input is a data structure of type SpectDG::Input containing the user defined parameters
{
	// see definitions in parameters.cpp
	EPS = input.EPS_ii ;
	JMAX = input.JMAX_ii ;
	JMIN = input.JMIN_ii ;
	pi = input.pi_ii ;
	omega = input.omega_ii ;
	rho = input.rho_ii ;
	f1 = input.f1_ii ;
	t1 = input.t1_ii ;
	t2 = input.t2_ii ;
	h = input.h_ii ;
	r1 = input.r1_ii ;
	M = input.M_ii ;
	N = input.N_ii ;
	M1 = input.M1_ii ;
	N1 = input.N1_ii ;
	dt0 = input.dt0_ii ;
	day = input.day_ii ;
	ah = input.ah_ii ;
	av = input.av_ii ;

	ncoef1=N*M*(4+N*M*(6+N*M));
	ncoef2=N1*M1*(4+N1*M1*(6+N1*M1));
	A=zeros<mat>(M*N+N1*M1,M*N+N1*M1);
	C=zeros<mat>(M*N+N1*M1,M*N+N1*M1);
	D=zeros<mat>(M*N+N1*M1,M*N+N1*M1);
	H=zeros<mat>(M*N+N1*M1,M*N+N1*M1);
	J=zeros<mat>(M*N+N1*M1,M*N+N1*M1);
	JJ=zeros<mat>(M*N+N1*M1,M*N+N1*M1);
	BW=zeros(M*N+N1*M1);
	BN=zeros<mat>(M*N+N1*M1);
	BNC=zeros<cube>(M*N+N1*M1,M*N+N1*M1,M*N+N1*M1);
	G=zeros<mat>(M*N+N1*M1,M*N+N1*M1);
	GG=zeros<mat>(M*N+N1*M1,M*N+N1*M1);
	h1=zeros(M*N+N1*M1); h2=zeros(M*N+N1*M1); h3=zeros(M*N+N1*M1);

	dt=dt0*24*60*60;//time step in seconds
	n=day*24*60*60/dt;//total number of time steps
	f2 = -f1;//final phi
	deltat=t2-t1;//theta range of the domain
	deltaf=f2-f1;//phi range of the domain
	l=2*pi/deltat;//spectral frequecny in theta direction
	d=pi/deltaf;//spectral frequecny in phi direction
	r2=r1+h;//earth radius at top of the ocean layer
	r_avg=0.5*r1+0.5*r2;//earth's average radius
	deltar=r2-r1;//ocean's total depth
	wa=0.9; wt=0.2;//wind asymmetry and tilt parameters
	t1_star=0.5*(t1+t2)-0.177*deltat-0.0647;//initial theta in small-scale
	t2_star=0.5*(t1+t2)+0.177*deltat-0.0647;//final theta in small-scale
	deltat_star=t2_star-t1_star;//theta range of small-scale box
	l_star=2*pi/deltat_star;//spectral frequecny in theta direction in small-scale
		
	vec Ydot1=zeros(M*N+N1*M1);vec Ydot2=zeros(M*N+N1*M1);vec Ydot3=zeros(M*N+N1*M1);vec Ydot4=zeros(M*N+N1*M1);
	vec Ydot=zeros(M*N+N1*M1);
	vec YYY=zeros(M*N+N1*M1);
	mat E=zeros<mat>(n);
	mat Y=zeros<mat>(M*N+N1*M1,n);
	t=0.0; f=0.0;//theta and phi coordinates
	FS1=0.0; FS2=0.0; FS3=0.0; FS4=0.0; FS5=0.0; FS6=0.0; FS7=0.0; FS8=0.0;
	H1=0.0; H2=0.0; H3=0.0;
	FN=0.0;
	Cs =0.1;//Smagorisnky constant
	nu_sgs1 = +pow(2.0, 0.5)*Cs*Cs*pow(deltaf*abs(cos(t1) - cos(t2))*pow(r_avg, 2)*h*1.0/30.0, 2.0 / 3.0);//Smagorinsky horizontal eddy viscosity
	nu_sgs2 = -pow(2.0, 0.5)*Cs*Cs*pow(deltaf*abs(cos(t1) - cos(t2))*pow(r_avg, 2)*h*1.0/30.0, 2.0 / 3.0);//Smagorinsky vertical eddy viscosity
	tau =- 2.0 * pi*0.8 / (10000.0 * 250.0);//wind stress amplitude
	
	tf=n*dt;
	y0 = zeros(1,M*N+N1*M1);

}