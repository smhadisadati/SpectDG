#include <armadillo>
#include <math.h>
#include <iostream>
#include <fstream>
#include <gnuplot_i.hpp>
#include <windows.h>
#include "globalvars.h"

using namespace std;
using namespace arma;

#define FUNC(x) ((*func)(x))

static double xsav, y1sav, y2sav;
static vec (*nrfunc)(double,double);
vec quad2d(vec (*func)(double, double), double, double, double, double);

vec coefFun1(double, double);
vec coefFun2(double, double);
colvec odeSys(double , colvec);
mat rksim(colvec (*func)(double, colvec), rowvec, double, double, double);

void spatialMat(void);
void yplot(mat);

double vphi1(double, double);
double vtheta1(double, double);
double weight_fun(int, int , double, double,double, double, double, double);
double time_der(int, int , double, double,double, double, double, double);
double visc_ver(int, int , double, double,double, double, double, double, double);
double visc_hor(int, int , double, double,double, double, double, double);
double eddy_ver(int, int , double, double,double, double, double, double, double);
double eddy_hor(int, int , double, double,double, double, double, double, double);
double Coriolis(int, int , double, double,double, double, double, double);


vec qtrapsim(vec (*func)(double), double a, double b, int num)
{
	vec s;
	for (int j=1;j<=JMAX;j++) {
		double x,tnm,del;
	vec s, sum;
	static vec s1, s2;
	num==1 ? s=s1 : s=s2;
	int it;
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
		
	}
	return s;
}