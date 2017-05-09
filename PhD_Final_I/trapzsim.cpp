#include "spectdg.h"

vec SpectDG::trapzsim(vec (SpectDG::*func)(double), double a, double b, int nt)
// simple fixed step trapazoidal integration for a faster solution
{
	double x,tnm,del;
	vec s, sum;
	int j;
	tnm=pow(2.0,nt-1);
	del=(b-a)/tnm;
	x=a+del;
	sum=0.5*(FUNC(a)+FUNC(b));
	for (j=2;j<=tnm;j++,x+=del) {
		sum = sum + FUNC(x);
	}
	s=sum*del;
	return s;
}