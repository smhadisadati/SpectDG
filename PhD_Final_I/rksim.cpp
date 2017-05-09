#include "spectdg.h"

mat SpectDG::rksim(colvec (SpectDG::*func)(double, colvec), rowvec y0, double t0, double tf, double th)
//fixed step forth order Runge-Kutta ODE solver in time
//inputs are the states derivative vector in the governing ODE system, vector of initial conditions, initial time, final time and time step respectively
//returns vector of states in the new time step
{
	colvec k1, k2, k3, k4, ys=y0.st();
	double ns=(tf-t0)/th+1;
	mat yout=zeros<mat>(ns,M*N+N1*M1);

	int i=0;
	for(double ti=t0;ti<=tf;ti+=th)
	{
		yout(i,span::all)=ys.st();

		k1=th*((this->*func)(ti,ys));
		k2=th*((this->*func)(ti+th/2,ys+k1/2));
		k3=th*((this->*func)(ti+th/2,ys+k2/2));
		k4=th*((this->*func)(ti+th,ys+k3));
		ys=ys+k1/6+k2/3+k3/3+k4/6;

		i++;
		//if (i%500==0) { // print out for debug
		//	cout<<i<<"\t"<<th<<"\t"<<ti<<endl<<"ys:"<<endl;
		//	ys.print();
		//}
	}
	return yout;
}