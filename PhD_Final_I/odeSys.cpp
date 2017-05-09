#include "spectdg.h"

colvec SpectDG::odeSys(double t , colvec Y)
//the state space form of the governing equations is built at time equal to t based on previously calculated spatial matrices and state vector at t
//inputs are the time and the state vector at that specific point of time
//returns the states derivative vector from the governing ODE system
{
	int k, i, j;
	vec YYY=Y(span(0,M*N+M1*N1-1),0);
	vec Ydot=zeros(M*N+M1*N1,1),
	
	BN=zeros(M*N+M1*N1);
	for(k=1;k<=M*N+M1*N1;k++){
		for(i=1;i<=M*N+M1*N1;i++){
			for(j=1;j<=M*N+N1*M1;j++){	
				BN(k-1)=BN(k-1)+BNC(i-1,j-1,k-1)*YYY(i-1)*YYY(j-1);	
			}
		}
	}

	for(i=1;i<=M*N+M1*N1;i++){
		for(j=1;j<=M*N+N1*M1;j++){
			G(i-1,j-1)=J(i-1,j-1)*abs(YYY(j-1));
			GG(i-1,j-1)=JJ(i-1,j-1)*abs(YYY(j-1));	
		}
	}
		
	vec B=BW*(250.0/deltar)+BN/r_avg+C*YYY/pow(r_avg,2.0)+D*YYY+H*YYY*(4000/deltar)+nu_sgs1*G*YYY/pow(r_avg,2.0)+nu_sgs2*GG*YYY/pow(h,2.0);
	
	
	Ydot=solve(A,B);

	

	return Ydot;
}