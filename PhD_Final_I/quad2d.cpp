#include "spectdg.h"

vec SpectDG::quad2d(vec (SpectDG::*func)(double, double), double x1, double x2, double y1, double y2)
//Returns the integral of a user-supplied function func over a two-dimensional region specified
//by the limits x1, x2, yy1 and yy2. Integration is performed by calling qtrap recursively.
//The necessary user-supplied functions have the following prototypes:
//double func(double x,double y,double z);
{
	y1sav = y1; y2sav = y2;
	nrfunc=func;
	return qtrap(&SpectDG::f11,x1,x2,1); // adaptive integration
	//return trapzsim(&SpectDG::f11,x1,x2,JMIN); // fixed step integration
}

vec SpectDG::f11(double xx)
{
	xsav=xx;
	return qtrap(&SpectDG::f22,y1sav,y2sav,2); // adaptive integration
	//return trapzsim(&SpectDG::f22,y1sav,y2sav,JMIN); // fixed step integration
}

vec SpectDG::f22(double yy) //The integrand f (x, y) evaluated at fixed x.
{
	return (this->*nrfunc)(xsav,yy);
}