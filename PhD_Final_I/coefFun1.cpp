#include "spectdg.h"

vec SpectDG::coefFun1(double t, double f)
//different terms of the governing equations are weighted and computed at a specific point of the large-scale solution domain
//inputs are theta and phi coordinates
//returns weighted terms of the governing equations at specific coordinates to be used in integrations and forming the spatial matirces
{
	vec coef = zeros(ncoef1);
	int ii=0, i, j, a, b, x, y;
	double V1, V2, V3, V4, V5, V6, FN, Z, tau_phi;
		
	for(j=1;j<=2*N-1;j+=2){
		for(i=1;i<=M;i++){
			V1=((-sin(j*d*(f-f1))*(cos(t)*sin(i*l*(t-t1))+i*l*cos(i*l*(t-t1))*sin(t))/(j*d))*cos(t)*cos(f)-(sin(i*l*(t-t1))*cos(j*d*(f-f1)))*sin(f))*sin(t);
			V2=((-sin(j*d*(f-f1))*(cos(t)*sin(i*l*(t-t1))+i*l*cos(i*l*(t-t1))*sin(t))/(j*d))*cos(t)*sin(f)+(sin(i*l*(t-t1))*cos(j*d*(f-f1)))*cos(f))*sin(t);
			V3=-(-sin(j*d*(f-f1))*(cos(t)*sin(i*l*(t-t1))+i*l*cos(i*l*(t-t1))*sin(t))/(j*d))*sin(t)*sin(t);
			if (t>t1_star & t<t2_star)
					{if ((t2-t)/deltat<(0.5+wt*((f-f1)/deltaf-0.5)))
						V4=0.0;//-sin(t)*(tau*wa)*sin((pi*(t2-t)/deltat)/(0.5+wt*((f-f1)/deltaf-0.5)))*weight_fun(i,j,t,f,t1,f1,l,d);//0.0
					else if ((t2-t)/deltat>=(0.5+wt*((f-f1)/deltaf-0.5)))
						V4=0.0;//sin(t)*(tau/wa)*sin((pi*((t2-t)/deltat-(0.5+wt*((f-f1)/deltaf-0.5))))/(0.5-wt*((f-f1)/deltaf-0.5)))*weight_fun(i,j,t,f,t1,f1,l,d);//0.0
					}
				else
					{if ((t2-t)/deltat<(0.5+wt*((f-f1)/deltaf-0.5)))
						V4=weight_fun(i,j,t,f,t1,f1,l,d)*(-sin(t)*(tau*wa)*sin((pi*(t2-t)/deltat)/(0.5+wt*((f-f1)/deltaf-0.5))));
					else if ((t2-t)/deltat>=(0.5+wt*((f-f1)/deltaf-0.5)))
						V4=weight_fun(i,j,t,f,t1,f1,l,d)*(sin(t)*(tau/wa)*sin((pi*((t2-t)/deltat-(0.5+wt*((f-f1)/deltaf-0.5))))/(0.5-wt*((f-f1)/deltaf-0.5))));
					}
			
			coef(ii) =pow(r_avg, 3.0)*deltar*V1; ii++;
			coef(ii) = pow(r_avg, 3)*deltar*V2; ii++;
			coef(ii) =pow(r_avg, 3)*deltar*V3; ii++;
			coef(ii)=V4; ii++;
			
			
			for(b=1;b<=2*N-1;b+=2){
				for(a=1;a<=M;a++){
					V1=time_der(a,b,t,f,t1,f1,l,d)*weight_fun(i,j,t,f,t1,f1,l,d);
					V2=visc_ver(a,b,t,f,t1,f1,l,d,r_avg)*weight_fun(i,j,t,f,t1,f1,l,d);
					V3=eddy_ver(a,b,t,f,t1,f1,l,d,r_avg)*visc_ver(a,b,t,f,t1,f1,l,d,r_avg)*weight_fun(i,j,t,f,t1,f1,l,d)/(av+2.0*av*h/r_avg);
					V4=visc_hor(a,b,t,f,t1,f1,l,d)*weight_fun(i,j,t,f,t1,f1,l,d);
					V5=eddy_hor(a,b,t,f,t1,f1,l,d,r_avg)*visc_hor(a,b,t,f,t1,f1,l,d)*weight_fun(i,j,t,f,t1,f1,l,d)/ah;
					V6=Coriolis(a,b,t,f,t1,f1,l,d)*weight_fun(i,j,t,f,t1,f1,l,d);
					coef(ii) = (V1); ii++;
					coef(ii) = (V2); ii++;
					coef(ii) = (V3); ii++;
					coef(ii) = (V4); ii++;
					coef(ii) = (V5); ii++;
					coef(ii) = (V6); ii++;
				}
			}
		}	
	}

	for(y=1;y<=2*N-1;y+=2){
		for(x=1;x<=M;x++){
			for(b=1;b<=2*N-1;b+=2){
				for(a=1;a<=M;a++){
					for(j=1;j<=2*N-1;j+=2){
						for(i=1;i<=M;i++){
							FN=sin(t)*(-((cos(b*d*(f-f1))*sin(d*j*(f-f1))*sin(a*l*(t-t1))*(sin(i*l*(t-t1))*cos(t)+i*l*cos(i*l*(t-t1))*sin(t)))/(d*j)+(cos(b*d*(f-f1))*sin(d*j*(f-f1))*sin(a*l*(t-t1))*(sin(i*l*(t-t1))*cos(t)+(i*i)*(l*l)*sin(i*l*(t-t1))*cos(t)*3.0+(i*i*i)*(l*l*l)*cos(i*l*(t-t1))*sin(t)+i*l*cos(i*l*(t-t1))*sin(t)*3.0))/(d*j)-b*d*cos(d*j*(f-f1))*sin(b*d*(f-f1))*sin(a*l*(t-t1))*1.0/pow(sin(t),2.0)*(sin(i*l*(t-t1))*cos(t)+i*l*cos(i*l*(t-t1))*sin(t))+(cos(b*d*(f-f1))*sin(d*j*(f-f1))*sin(a*l*(t-t1))*1.0/tan(t)*(sin(i*l*(t-t1))*sin(t)+(i*i)*(l*l)*sin(i*l*(t-t1))*sin(t)-i*l*cos(i*l*(t-t1))*cos(t)*2.0)*2.0)/(d*j)-((b*b)*d*cos(b*d*(f-f1))*sin(d*j*(f-f1))*sin(a*l*(t-t1))*1.0/pow(sin(t),2.0)*(sin(i*l*(t-t1))*cos(t)+i*l*cos(i*l*(t-t1))*sin(t)))/j+(a*l*cos(b*d*(f-f1))*sin(d*j*(f-f1))*cos(a*l*(t-t1))*(sin(i*l*(t-t1))*sin(t)+(i*i)*(l*l)*sin(i*l*(t-t1))*sin(t)-i*l*cos(i*l*(t-t1))*cos(t)*2.0))/(d*j)-(a*l*cos(b*d*(f-f1))*sin(d*j*(f-f1))*cos(a*l*(t-t1))*1.0/tan(t)*(sin(i*l*(t-t1))*cos(t)+i*l*cos(i*l*(t-t1))*sin(t)))/(d*j))-((a*d*j*l*cos(b*d*(f-f1))*sin(d*j*(f-f1))*cos(a*l*(t-t1))*sin(i*l*(t-t1)))/sin(t)+(d*i*j*l*cos(b*d*(f-f1))*sin(d*j*(f-f1))*cos(i*l*(t-t1))*sin(a*l*(t-t1)))/sin(t))+((cos(d*j*(f-f1))*sin(b*d*(f-f1))*(sin(i*l*(t-t1))*cos(t)+i*l*cos(i*l*(t-t1))*sin(t))*(sin(a*l*(t-t1))*sin(t)+(a*a)*(l*l)*sin(a*l*(t-t1))*sin(t)-a*l*cos(a*l*(t-t1))*cos(t)*2.0))/(b*d*sin(t))+(cos(d*j*(f-f1))*sin(b*d*(f-f1))*(sin(a*l*(t-t1))*cos(t)+a*l*cos(a*l*(t-t1))*sin(t))*(sin(i*l*(t-t1))*sin(t)+(i*i)*(l*l)*sin(i*l*(t-t1))*sin(t)-i*l*cos(i*l*(t-t1))*cos(t)*2.0))/(b*d*sin(t))-(cos(d*j*(f-f1))*sin(b*d*(f-f1))*1.0/tan(t)*(sin(a*l*(t-t1))*cos(t)+a*l*cos(a*l*(t-t1))*sin(t))*(sin(i*l*(t-t1))*cos(t)+i*l*cos(i*l*(t-t1))*sin(t))*2.0)/(b*d*sin(t))));
							Z=FN*(d*y*sin(d*y*(f-f1))*sin(l*x*(t-t1))-(sin(d*y*(f-f1))*cos(t)*(sin(l*x*(t-t1))*cos(t)+l*x*cos(l*x*(t-t1))*sin(t)))/(d*y)+(sin(d*y*(f-f1))*sin(t)*(sin(l*x*(t-t1))*sin(t)+(l*l)*(x*x)*sin(l*x*(t-t1))*sin(t)-l*x*cos(l*x*(t-t1))*cos(t)*2.0))/(d*y))/sin(t);
							coef(ii)=(Z); ii++;
						}
					}
				}
			}
		}
	}
		
	return coef;
}