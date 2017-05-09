=========================================================
Credits:
=========================================================
Authors: S. Elnaz Naghibi(1)* and S.M.Hadi Sadati(2)
(1)*s.e.naghibi@qmul.ac.uk
School of Engineering and Materials Science, Queen Mary University of London
(2)seyedmohammadhadi.sadati@kcl.ac.uk
Informatics Department, King's College London
Last update: 21/08/2016

=========================================================
General Notes and Installation:
=========================================================
This package presents a two-scale spectral solution for the double-gyre problem in a  
turbulent regime.
The codes are written and tested using Microsoft Visual C++ 2010 Express.
armadillo, GNUplot, Boost, BLAS and LAPACK libraries for C++ are necessary to be installed for this header file content to work properly.
Please use /MP /O2 /Zi /Gm- command-line options for a faster and compatible run-time.

=========================================================
Files description:
=========================================================
>> spectdg.h
Header file for global definitions, "SpectDG" class and declaration of "parameters"  
function

>> run.cpp
Spectral solution procedure
This method runs the space and time integrations.

>> parameters.cpp
Spectral solution and flow parameters are assigned by the user.
functions:
SpectDG::Input parameters
It creates the input data structure based on user defined parameters and returns a data structure of type SpectDG::Input.

>> initialize.cpp
It sets the input parameters for a SpectDG object and initializes the solution parameters.
Input is a data structure of type SpectDG::Input containing the user defined parameters.

>> main.cpp
It gathers user input parameters, instantiates a spectral solution object, runs the spectral solution and saves and plots the solution results by calling the corresponding functions and methods.

>> governing_equations_terms.cpp
Inputs of the following functions are the double summation indices, theta and phi  
coordinates, initial theta and phi, meridional and zonal frequencies
of the spectral solutions and average radius of the ocean layer. Outputs are weighting  
functions, time derivatives, vertical and horizontal viscosity terms, vertical and  
horizontal eddy viscosity terms, and Coriolis terms.
methods:
SpectDG::weight_fun, SpectDG::time_der, SpectDG::visc_ver, SpectDG::visc_hor, 
SpectDG::eddy_ver, SpectDG::eddy_hor, SpectDG::Coriolis

>> coefFun1.cpp
Different terms of the governing equations are weighted and computed at a specific point 
of the large-scale solution domain.
Inputs are theta and phi coordinates.
It returns weighted terms of the governing equations at specific coordinates to be used in integrations and forming the spatial matrices.

>> coefFun2.cpp
Different terms of the governing equations are weighted and computed at a specific point 
of the small-scale solution domain.
Inputs are theta and phi coordinates.
It returns weighted terms of the governing equations at specific coordinates to be used in integrations and forming the spatial matrices.

>> trapzsim.cpp
 Simple fixed step trapezoidal integration for a faster solution

>> quad2d.cpp
It returns the integral of a user-supplied function func over a two-dimensional region  
specified by the limits x1, x2, yy1 and yy2. Integration is performed by calling qtrap  
recursively.
methods:
SpectDG::f11, SpectDG::f22

>> qtrap.cpp
It returns the integral of the function func from a to b. The parameters EPS can be set to the desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum  
allowed number of steps. Integration is performed by the trapezoidal rule.
methods:
SpectDG::trapzd
This routine computes the nth stage of refinement of an extended trapezoidal rule. func is input as a pointer to the function to be integrated between limits a and b, also input. When called with n=1, the routine returns the crudest estimate of a b f (x)dx. Subsequent calls with n=2,3,... (in that sequential order) will improve the accuracy by adding 2n-2 additional interior points.

>> spatialMat.cpp
System ODE matrices are built by integrating weighted parts of the governing equations  
computed in all points.

>> odeSys.cpp
The state space form of the governing equations is built at time equal to t based on  
previously calculated spatial matrices and state vector at t.
Inputs are the time and the state vector at that specific point of time.
It returns the states derivative vector from the governing ODE system.

>> rksim.cpp
Fixed step forth order Runge-Kutta ODE solver in time
Inputs are the states derivative vector in the governing ODE system, vector of initial  
conditions, initial time, final time and time step respectively.
It returns vector of states in the new time step.

>> output.cpp
Output method saves and plots the time solution of the simulation.

>> yplot.cpp
It plots the first velocity amplitude in large-scale and small-scale vs. time step.
