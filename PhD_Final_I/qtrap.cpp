#include "spectdg.h"

vec SpectDG::qtrap(vec (SpectDG::*func)(double), double a, double b, int num)
//Returns the integral of the function func from a to b. The parameters EPS can be set to the
//desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum allowed
//number of steps. Integration is performed by the trapezoidal rule.
{
	int j;
	vec s, olds=zeros(ncoef); //Initial value of olds is arbitrary.
	boolean err;
	for (j=1;j<=JMAX;j++) {
		s=trapzd(func,a,b,j,num);
		if (j > JMIN) //Avoid spurious early convergence.
		{
			//cout<<j<<endl<<endl;
			err=true;
			for(int i=0;i<ncoef;i++) if(abs(s(i)-olds(i))/abs(olds(i))>EPS) {
				//cout << i << "\t" << abs(s(i)-olds(i))/abs(olds(i)) << endl;
				err=false;
			}
			if (err) {
				cout<<j<<endl;
				return s;
			}
		}
		olds=s;
	}
	if (JMAX != JMIN) cout << "Too many steps in routine qtrap" << endl;
	return s;
}

vec SpectDG::trapzd(vec (SpectDG::*func)(double), double a, double b, int nt, int num)
//This routine computes the nth stage of refinement of an extended trapezoidal rule. func is input
//as a pointer to the function to be integrated between limits a and b, also input. When called with
//n=1, the routine returns the crudest estimate of a b f (x)dx. Subsequent calls with n=2,3,...
//(in that sequential order) will improve the accuracy by adding 2n-2 additional interior points.
{
	double x,tnm,del;
	vec s, sum;
	static vec s1, s2;
	num==1 ? s=s1 : s=s2;
	int it,j;
	if (nt == 1) {
		s=0.5*(b-a)*(FUNC(a)+FUNC(b));
	} else {
		for (it=1,j=1;j<nt-1;j++)
			it <<= 1;
		tnm=it;
		del=(b-a)/tnm; //This is the spacing of the points to be added.
		x=a+0.5*del;
		for (sum=zeros(ncoef),j=1;j<=it;j++,x+=del) sum = sum + FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm); //This replaces s by its refined value.
	}
	num==1 ? s1=s : s2=s;
	return s;
}