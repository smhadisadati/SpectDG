#include "spectdg.h"
/*Inputs of the following functions are the double summation indices, theta and phi coordinates, initial theta and phi, meridional and zonal frequencies
of the spectral solutions and average radius of the ocean layer*/

double SpectDG::weight_fun(int i, int j, double t, double f,double t_1, double f_1, double l,double d)
//weighting functions
{
	double wf=((d*j*sin(d*j*(f-f_1))*sin(i*l*(t-t_1))-(sin(d*j*(f-f_1))*cos(t)*(sin(i*l*(t-t_1))*cos(t)+i*l*cos(i*l*(t-t_1))*sin(t)))/(d*j)+(sin(d*j*(f-f_1))*sin(t)*(sin(i*l*(t-t_1))*sin(t)+(i*i)*(l*l)*sin(i*l*(t-t_1))*sin(t)-i*l*cos(i*l*(t-t_1))*cos(t)*2.0))/(d*j))/sin(t));
	return wf;
}

double SpectDG::time_der(int i, int j, double t, double f,double t_1, double f_1, double l,double d)
//time derivative terms
{
	double tder=((d*j*sin(d*j*(f-f_1))*sin(i*l*(t-t_1))-(sin(d*j*(f-f_1))*cos(t)*(sin(i*l*(t-t_1))*cos(t)+i*l*cos(i*l*(t-t_1))*sin(t)))/(d*j)+(sin(d*j*(f-f_1))*sin(t)*(sin(i*l*(t-t_1))*sin(t)+(i*i)*(l*l)*sin(i*l*(t-t_1))*sin(t)-i*l*cos(i*l*(t-t_1))*cos(t)*2.0))/(d*j))/sin(t))*sin(t);
	return tder;
}

double SpectDG::visc_ver(int i, int j, double t, double f,double t_1, double f_1, double l,double d, double r_avg)
//vertical viscosity terms
{
	double viscv=(av+2.0*av*h/r_avg)*((d*j*sin(d*j*(f-f_1))*sin(i*l*(t-t_1))-(sin(d*j*(f-f_1))*cos(t)*(sin(i*l*(t-t_1))*cos(t)+i*l*cos(i*l*(t-t_1))*sin(t)))/(d*j)+(sin(d*j*(f-f_1))*sin(t)*(sin(i*l*(t-t_1))*sin(t)+(i*i)*(l*l)*sin(i*l*(t-t_1))*sin(t)-i*l*cos(i*l*(t-t_1))*cos(t)*2.0))/(d*j))/sin(t))*sin(t);
	return viscv;
}



double SpectDG::visc_hor(int a, int b, double t, double f,double t_1, double f_1, double l,double d)
//horizontal viscosity terms
{
	double FS1 =-(b*d*sin(b*d*(f-f_1))*sin(a*l*(t-t_1))*sin(t)+(a*a)*b*d*(l*l)*sin(b*d*(f-f_1))*sin(a*l*(t-t_1))*sin(t)-a*b*d*l*sin(b*d*(f-f_1))*cos(a*l*(t-t_1))*cos(t)*2.0)/sin(t)-cos(t)*1.0/pow(sin(t),2.0)*(b*d*sin(b*d*(f-f_1))*sin(a*l*(t-t_1))*cos(t)+a*b*d*l*sin(b*d*(f-f_1))*cos(a*l*(t-t_1))*sin(t));
	double FS2 =-(b*b*b)*(d*d*d)*sin(b*d*(f-f_1))*sin(a*l*(t-t_1))*1.0/pow(sin(t),2.0);
	double FS3 =(b*d*sin(b*d*(f-f_1))*1.0/tan(t)*(sin(a*l*(t-t_1))*cos(t)+a*l*cos(a*l*(t-t_1))*sin(t)))/sin(t);
	double FS4 = -sin(t)*(((sin(b*d*(f-f_1))*cos(t)*(sin(a*l*(t-t_1))*cos(t)+a*l*cos(a*l*(t-t_1))*sin(t)))/(b*d)-(sin(b*d*(f-f_1))*sin(t)*(sin(a*l*(t-t_1))*sin(t)+(a*a)*(l*l)*sin(a*l*(t-t_1))*sin(t)-a*l*cos(a*l*(t-t_1))*cos(t)*2.0))/(b*d))/sin(t)-((sin(b*d*(f-f_1))*cos(t)*(sin(a*l*(t-t_1))*cos(t)+a*l*cos(a*l*(t-t_1))*sin(t)))/(b*d)+(sin(b*d*(f-f_1))*cos(t)*(sin(a*l*(t-t_1))*cos(t)+(a*a)*(l*l)*sin(a*l*(t-t_1))*cos(t)*3.0+(a*a*a)*(l*l*l)*cos(a*l*(t-t_1))*sin(t)+a*l*cos(a*l*(t-t_1))*sin(t)*3.0)*3.0)/(b*d)-(sin(b*d*(f-f_1))*sin(t)*(sin(a*l*(t-t_1))*sin(t)+(a*a)*(l*l)*sin(a*l*(t-t_1))*sin(t)-a*l*cos(a*l*(t-t_1))*cos(t)*2.0)*3.0)/(b*d)-(sin(b*d*(f-f_1))*sin(t)*(sin(a*l*(t-t_1))*sin(t)+(a*a)*(l*l)*sin(a*l*(t-t_1))*sin(t)*6.0+(a*a*a*a)*(l*l*l*l)*sin(a*l*(t-t_1))*sin(t)-a*l*cos(a*l*(t-t_1))*cos(t)*4.0-(a*a*a)*(l*l*l)*cos(a*l*(t-t_1))*cos(t)*4.0))/(b*d))/sin(t)+cos(t)*1.0/pow(sin(t),2.0)*((sin(b*d*(f-f_1))*sin(t)*(sin(a*l*(t-t_1))*cos(t)+a*l*cos(a*l*(t-t_1))*sin(t)))/(b*d)+(sin(b*d*(f-f_1))*cos(t)*(sin(a*l*(t-t_1))*sin(t)+(a*a)*(l*l)*sin(a*l*(t-t_1))*sin(t)-a*l*cos(a*l*(t-t_1))*cos(t)*2.0)*2.0)/(b*d)+(sin(b*d*(f-f_1))*sin(t)*(sin(a*l*(t-t_1))*cos(t)+(a*a)*(l*l)*sin(a*l*(t-t_1))*cos(t)*3.0+(a*a*a)*(l*l*l)*cos(a*l*(t-t_1))*sin(t)+a*l*cos(a*l*(t-t_1))*sin(t)*3.0))/(b*d))*2.0+pow(cos(t),2.0)*1.0/pow(sin(t),3.0)*((sin(b*d*(f-f_1))*cos(t)*(sin(a*l*(t-t_1))*cos(t)+a*l*cos(a*l*(t-t_1))*sin(t)))/(b*d)-(sin(b*d*(f-f_1))*sin(t)*(sin(a*l*(t-t_1))*sin(t)+(a*a)*(l*l)*sin(a*l*(t-t_1))*sin(t)-a*l*cos(a*l*(t-t_1))*cos(t)*2.0))/(b*d))*2.0);
	double FS5 =-(b*d*sin(b*d*(f-f_1))*(sin(a*l*(t-t_1))*sin(t)+(a*a)*(l*l)*sin(a*l*(t-t_1))*sin(t)-a*l*cos(a*l*(t-t_1))*cos(t)*2.0))/sin(t);
	double FS6 =a*b*d*l*sin(b*d*(f-f_1))*cos(a*l*(t-t_1))*1.0/tan(t)*-2.0;
	double FS7=cos(t)*(((sin(b*d*(f-f_1))*sin(t)*(sin(a*l*(t-t_1))*cos(t)+a*l*cos(a*l*(t-t_1))*sin(t)))/(b*d)+(sin(b*d*(f-f_1))*cos(t)*(sin(a*l*(t-t_1))*sin(t)+(a*a)*(l*l)*sin(a*l*(t-t_1))*sin(t)-a*l*cos(a*l*(t-t_1))*cos(t)*2.0)*2.0)/(b*d)+(sin(b*d*(f-f_1))*sin(t)*(sin(a*l*(t-t_1))*cos(t)+(a*a)*(l*l)*sin(a*l*(t-t_1))*cos(t)*3.0+(a*a*a)*(l*l*l)*cos(a*l*(t-t_1))*sin(t)+a*l*cos(a*l*(t-t_1))*sin(t)*3.0))/(b*d))/sin(t)+cos(t)*1.0/pow(sin(t),2.0)*((sin(b*d*(f-f_1))*cos(t)*(sin(a*l*(t-t_1))*cos(t)+a*l*cos(a*l*(t-t_1))*sin(t)))/(b*d)-(sin(b*d*(f-f_1))*sin(t)*(sin(a*l*(t-t_1))*sin(t)+(a*a)*(l*l)*sin(a*l*(t-t_1))*sin(t)-a*l*cos(a*l*(t-t_1))*cos(t)*2.0))/(b*d)));
	double FS8=b*d*sin(b*d*(f-f_1))*sin(a*l*(t-t_1))*1.0/pow(sin(t),2.0)*2.0;

	double visch=ah*(FS1+FS2+FS3+FS4+FS5+FS6+FS7+FS8); 
	return visch;
}


double SpectDG::eddy_ver(int a, int b, double t, double f,double t_1, double f_1, double l,double d, double r_avg)
//vertical eddy viscosity terms
{
	
	double FS7=sin(a*l*(t-t_1))*cos(b*d*(f-f_1))*(1.0/h-1.0/r_avg);
	double FS8=-sin(b*d*(f-f_1))*(cos(t)*sin(a*l*(t-t_1))+a*l*sin(t)*cos(a*l*(t-t_1)))*(1.0/h-1.0/r_avg)/(b*d);
	double eddyv=0.5*pow(2*FS7*FS7+2*FS8*FS8,0.5);

	return eddyv;
}

double SpectDG::eddy_hor(int a, int b, double t, double f,double t_1, double f_1, double l,double d, double r_avg)
//horizontal eddy viscosity terms
{
	
	double FS1=a*l*cos(b*d*(f-f_1))*cos(a*l*(t-t_1));
	double FS2=-(b*d*sin(b*d*(f-f_1))*sin(a*l*(t-t_1)))/sin(t);
	double FS3=(sin(b*d*(f-f_1))*1.0/tan(t)*(sin(a*l*(t-t_1))*cos(t)+a*l*cos(a*l*(t-t_1))*sin(t)))/(b*d);
	double FS4=(sin(b*d*(f-f_1))*(sin(a*l*(t-t_1))*sin(t)+(a*a)*(l*l)*sin(a*l*(t-t_1))*sin(t)-a*l*cos(a*l*(t-t_1))*cos(t)*2.0))/(b*d);
	double FS5=-(cos(b*d*(f-f_1))*(sin(a*l*(t-t_1))*cos(t)+a*l*cos(a*l*(t-t_1))*sin(t)))/sin(t);
	double FS6=cos(b*d*(f-f_1))*sin(a*l*(t-t_1))*1.0/tan(t);

	double eddyh=0.5*pow((2.0*FS1*FS1+2.0*(FS2+FS3+FS4)*(FS2+FS3+FS4)+2.0*(FS5+FS6)*(FS5+FS6))/(r_avg*r_avg),0.5);

	return eddyh;
}

double SpectDG::Coriolis(int a, int b, double t, double f,double t_1, double f_1, double l,double d)
//Coriolis terms
{
	double coriolis=2.0*omega*(sin(a*l*(t-t_1))*cos(b*d*(f-f_1)))*sin(t)*sin(t);
	return coriolis;
}
