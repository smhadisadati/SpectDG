#include <gnuplot_i.hpp> // needs GNUplot to be installed and relevant system variables to be  set
#include "spectdg.h"

void SpectDG::yplot(mat yout)
//plots the first velocity amplitude in large-scale and small-scale vs. time step
{
	
	cout<<"Save and plot..."<<endl;
	vec Yout1 = yout(span::all,0); // vector data prepration 
	vec Yout2 = yout(span::all,M*N);

	static Gnuplot g1 ("Plot"); // plot object
	g1.set_style("lines").set_grid(); // set style
	g1.set_xlabel ("time steps") .set_ylabel ("Y");
	g1.plot_x ( Yout1 ,"Y1_Large_scale"); // plot
	g1.plot_x ( Yout2 ,"Y1_small_scale");

}
