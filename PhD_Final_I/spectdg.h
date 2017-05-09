/*
In The name of God

This package presents a two-scale spectral solution for the double-gyre problem in a  
turbulent regime. The codes are written and tested using Microsoft Visual C++ 2010 Express.
armadillo, GNUplot, Boost, BLAS and LAPACK libraries for C++ are necessary to be installed for this header file content to work properly.
Please use /MP /O2 /Zi /Gm- command-line options for a faster and compatible run-time.

Authors: S. Elnaz Naghibi(1)* and S.M.Hadi Sadati(2)
(1)*s.e.naghibi@qmul.ac.uk
School of Engineering and Materials Science, Queen Mary University of London
(2)seyedmohammadhadi.sadati@kcl.ac.uk
Informatics Department, King's College London
Last update: 21/08/2016

*/


// header files
#include <armadillo> // linear algebra package
#include <math.h>
#include <iostream>
#include <fstream>
#include <windows.h>

using namespace std;
using namespace arma; // linear algebra namespace

#ifndef GLOBALVARS_H
#define GLOBALVARS_H

// global definitions
#define FUNC(x) ((this->*func)(x))

// class definition
class SpectDG {
// SpectDG class for a spectral solution object
private:
	// variables (see definitions in parameters.cpp and initialize.cpp)
	double EPS ;
	int JMAX ;
	int JMIN ;

	double pi ;
	double omega ;
	double rho ;
	double f1 ;
	double t1 ;
	double t2 ;
	double h ;
	double r1 ;
	int M ;
	int N ;
	int M1 ;
	int N1 ;
	double dt0 ;

	int day ;
	double ah ;
	double av ;

	double xsav, y1sav, y2sav;
	int ncoef, ncoef1, ncoef2;
	mat A;
	mat C;
	mat D;
	mat H;
	mat J;
	mat JJ;
	vec BW;
	vec BN;
	cube BNC;
	mat G;
	mat GG;
	double L1, L2, L3;
	vec h1, h2, h3;

	double dt;
	int n;
	double f2;
	double deltat;
	double deltaf;
	double l;
	double d;
	double r2;
	double r_avg;
	double deltar;
	double wa, wt;
	double t1_star;
	double t2_star;
	double deltat_star;
	double l_star;
	
	vec Ydot1, Ydot2, Ydot3, Ydot4;
	vec Ydot;
	vec YYY;
	mat E;
	mat Y;
	double t, f;
	double FS1, FS2, FS3, FS4, FS5, FS6, FS7, FS8;
	double H1, H2, H3;
	double FN;
	double Cs;
	double nu_sgs1;
	double nu_sgs2;
	double tau;
	
	double tf;
	rowvec y0;
	mat yout;

	vec (SpectDG::*nrfunc)(double,double); // stors pointer to the integrand method

	// private function declaration
	// see relevant .cpp files for more information
	vec quad2d(vec (SpectDG::*)(double, double), double, double, double, double);

	vec qtrap(vec (SpectDG::*)(double), double a, double b, int num);
	vec f11(double xx);
	vec f22(double yy);
	vec trapzd(vec (SpectDG::*)(double), double a, double b, int n, int num);
	vec trapzsim(vec (SpectDG::*)(double), double, double, int);

	vec coefFun1(double, double);
	vec coefFun2(double, double);
	colvec odeSys(double , colvec);
	mat rksim(colvec (SpectDG::*func)(double, colvec), rowvec, double, double, double);

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
	
public:
	//input parameters structure
	struct Input { // stores the user input parameters
		double EPS_ii ;
		int JMAX_ii ;
		int JMIN_ii ;

		double pi_ii ;
		double omega_ii ;
		double rho_ii ;
		double f1_ii ;
		double t1_ii ;
		double t2_ii ;
		double h_ii ;
		double r1_ii ;
		int M_ii ;
		int N_ii ;
		int M1_ii ;
		int N1_ii ;
		double dt0_ii ;

		int day_ii ;
		double ah_ii ;
		double av_ii ;
	};

	//constructor:
	// creates an instance of the class and initiaize the object parameters
	SpectDG(SpectDG::Input input)
	{
		cout<<"Initialization..."<<endl;
		initialize(input); // initialization
	}

	// destructor
	~SpectDG(){}

	// public function declarations:
	void initialize(SpectDG::Input); // initialization
	void run(); // runs solution
	void output(string); // saves and plots output
	void printInput(); // for debug
};

// global function declaration
SpectDG::Input parameters(); // retrieves user input parameters

#endif