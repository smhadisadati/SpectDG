#include "spectdg.h"

void SpectDG::run()
// spectral solution procedure
// this method runs the space and time integrations
{
	cout<<"Spatial matrices integration..."<<endl;
	spatialMat(); // space integration
	cout<<"Time solution begins..."<<endl;
	yout = rksim(&SpectDG::odeSys, y0, 0.0, tf, dt); // time integration
}