#include "spectdg.h" // include the class header file

int main()
{
	// spectral solution procedure
	// Initialization
	SpectDG::Input input = parameters();  // gathers user input parameters
	SpectDG spectDG(input); // instantiates a spectral solution object
	//spectDG.printInput(); // for debug
	spectDG.run(); // runs the spectral solution
	spectDG.output("yout.TXT"); // saves and plots the solution results
		
	return 0;
}

