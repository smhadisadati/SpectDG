#include "spectdg.h"

void SpectDG::output(string name)
// output method saves and plots the time solution of the simulation
//input is a string for the name and location of the save file
{
	yout.save(name, arma_ascii); // save solution
	yplot(yout); // plot output
	
	cout<<"Done!"<<endl<<"Press enter to terminate."<<endl;
	getchar();
}

void SpectDG::printInput()
// a method to prints the input parameters for debug
{
	cout<<EPS<<"\t"<<JMAX<<"\t"<<JMIN<<"\t"<<pi<<"\t"<<omega<<"\t"<<rho<<"\t"<<f1<<"\t"<<t1<<"\t"<<t2<<"\t"<<h<<"\t"<<r1<<"\t"<<M<<"\t"<<N<<"\t"<<M1<<"\t"<<N1<<"\t"<<dt0<<"\t"<<day<<"\t"<<ah<<"\t"<<av<<endl;
}