#include "spectdg.h"

vec SpectDG::coefFun2(double t, double f)
//different terms of the governing equations are weighted and computed at a specific point of the small-scale solution domain
//inputs are theta and phi coordinates
//returns weighted terms of the governing equations at specific coordinates to be used in integrations and forming the spatial matirces
{
	vec coef = zeros(ncoef2);
	int ii=0, i, j, a, b, x, y;
	double V1, V2, V3, V4, V5, V6, FN, Z, tau_phi;
	
	for(j=1;j<=2*N1-1;j+=2){
		for(i=1;i<=M1;i++){
			V1=((-sin(j*d*(f-f1))*(cos(t)*sin(i*l_star*(t-t1_star))+i*l_star*cos(i*l_star*(t-t1_star))*sin(t))/(j*d))*cos(t)*cos(f)-(sin(i*l_star*(t-t1_star))*cos(j*d*(f-f1)))*sin(f))*sin(t);
			V2=((-sin(j*d*(f-f1))*(cos(t)*sin(i*l_star*(t-t1_star))+i*l_star*cos(i*l_star*(t-t1_star))*sin(t))/(j*d))*cos(t)*sin(f)+(sin(i*l_star*(t-t1_star))*cos(j*d*(f-f1)))*cos(f))*sin(t);
			V3=-(-sin(j*d*(f-f1))*(cos(t)*sin(i*l_star*(t-t1_star))+i*l_star*cos(i*l_star*(t-t1_star))*sin(t))/(j*d))*sin(t)*sin(t);
			if ((t2-t)/deltat<(0.5+wt*((f-f1)/deltaf-0.5)))
				V4=-1.0*sin(t)*(tau*wa)*sin((pi*(t2-t)/deltat)/(0.5+wt*((f-f1)/deltaf-0.5)))*weight_fun(i,j,t,f,t1_star,f1,l_star,d);
			else if ((t2-t)/deltat>=(0.5+wt*((f-f1)/deltaf-0.5)))
				V4=1.0*sin(t)*(tau/wa)*sin((pi*((t2-t)/deltat-(0.5+wt*((f-f1)/deltaf-0.5))))/(0.5-wt*((f-f1)/deltaf-0.5)))*weight_fun(i,j,t,f,t1_star,f1,l_star,d);					
			coef(ii) = pow(r_avg, 3)*deltar*(V1); ii++;
			coef(ii) = pow(r_avg, 3)*deltar*(V2); ii++;
			coef(ii) = pow(r_avg, 3)*deltar*(V3); ii++;
			coef(ii)= (V4); ii++;
			

			for(b=1;b<=2*N1-1;b+=2){
				for(a=1;a<=M1;a++){
					V1=time_der(a,b,t,f,t1_star,f1,l_star,d)*weight_fun(i,j,t,f,t1_star,f1,l_star,d);
					V2=visc_ver(a,b,t,f,t1_star,f1,l_star,d,r_avg)*weight_fun(i,j,t,f,t1_star,f1,l_star,d);
					V3=eddy_ver(a,b,t,f,t1_star,f1,l_star,d,r_avg)*visc_ver(a,b,t,f,t1_star,f1,l_star,d,r_avg)*weight_fun(i,j,t,f,t1_star,f1,l_star,d)/(av+2.0*av*h/r_avg);
					V4=visc_hor(a,b,t,f,t1_star,f1,l_star,d)*weight_fun(i,j,t,f,t1_star,f1,l_star,d);
					V5=eddy_hor(a,b,t,f,t1_star,f1,l_star,d,r_avg)*visc_hor(a,b,t,f,t1_star,f1,l_star,d)*weight_fun(i,j,t,f,t1_star,f1,l_star,d)/ah;
					V6=Coriolis(a,b,t,f,t1_star,f1,l_star,d)*weight_fun(i,j,t,f,t1_star,f1,l_star,d);
					coef(ii) =  (V1) ; ii++;
					coef(ii) = (V2) ; ii++;
					coef(ii) = (V3) ; ii++;
					coef(ii) = (V4) ; ii++;
					coef(ii) = (V5) ; ii++;
					coef(ii) = (V6) ; ii++;
				}
				
			}
		}
	}

	for(y=1;y<=2*N1-1;y+=2){
		for(x=1;x<=M1;x++){
			for(b=1;b<=2*N1-1;b+=2){
				for(a=1;a<=M1;a++){
					for(j=1;j<=2*N1-1;j+=2){
						for(i=1;i<=M1;i++){
							FN=sin(t)*(-((cos(b*d*(f-f1))*sin(d*j*(f-f1))*sin(a*l_star*(t-t1_star))*(sin(i*l_star*(t-t1_star))*cos(t)+i*l_star*cos(i*l_star*(t-t1_star))*sin(t)))/(d*j)+(cos(b*d*(f-f1))*sin(d*j*(f-f1))*sin(a*l_star*(t-t1_star))*(sin(i*l_star*(t-t1_star))*cos(t)+(i*i)*(l_star*l_star)*sin(i*l_star*(t-t1_star))*cos(t)*3.0+(i*i*i)*(l_star*l_star*l_star)*cos(i*l_star*(t-t1_star))*sin(t)+i*l_star*cos(i*l_star*(t-t1_star))*sin(t)*3.0))/(d*j)-b*d*cos(d*j*(f-f1))*sin(b*d*(f-f1))*sin(a*l_star*(t-t1_star))*1.0/pow(sin(t),2.0)*(sin(i*l_star*(t-t1_star))*cos(t)+i*l_star*cos(i*l_star*(t-t1_star))*sin(t))+(cos(b*d*(f-f1))*sin(d*j*(f-f1))*sin(a*l_star*(t-t1_star))*1.0/tan(t)*(sin(i*l_star*(t-t1_star))*sin(t)+(i*i)*(l_star*l_star)*sin(i*l_star*(t-t1_star))*sin(t)-i*l_star*cos(i*l_star*(t-t1_star))*cos(t)*2.0)*2.0)/(d*j)-((b*b)*d*cos(b*d*(f-f1))*sin(d*j*(f-f1))*sin(a*l_star*(t-t1_star))*1.0/pow(sin(t),2.0)*(sin(i*l_star*(t-t1_star))*cos(t)+i*l_star*cos(i*l_star*(t-t1_star))*sin(t)))/j+(a*l_star*cos(b*d*(f-f1))*sin(d*j*(f-f1))*cos(a*l_star*(t-t1_star))*(sin(i*l_star*(t-t1_star))*sin(t)+(i*i)*(l_star*l_star)*sin(i*l_star*(t-t1_star))*sin(t)-i*l_star*cos(i*l_star*(t-t1_star))*cos(t)*2.0))/(d*j)-(a*l_star*cos(b*d*(f-f1))*sin(d*j*(f-f1))*cos(a*l_star*(t-t1_star))*1.0/tan(t)*(sin(i*l_star*(t-t1_star))*cos(t)+i*l_star*cos(i*l_star*(t-t1_star))*sin(t)))/(d*j))-((a*d*j*l_star*cos(b*d*(f-f1))*sin(d*j*(f-f1))*cos(a*l_star*(t-t1_star))*sin(i*l_star*(t-t1_star)))/sin(t)+(d*i*j*l_star*cos(b*d*(f-f1))*sin(d*j*(f-f1))*cos(i*l_star*(t-t1_star))*sin(a*l_star*(t-t1_star)))/sin(t))+((cos(d*j*(f-f1))*sin(b*d*(f-f1))*(sin(i*l_star*(t-t1_star))*cos(t)+i*l_star*cos(i*l_star*(t-t1_star))*sin(t))*(sin(a*l_star*(t-t1_star))*sin(t)+(a*a)*(l_star*l_star)*sin(a*l_star*(t-t1_star))*sin(t)-a*l_star*cos(a*l_star*(t-t1_star))*cos(t)*2.0))/(b*d*sin(t))+(cos(d*j*(f-f1))*sin(b*d*(f-f1))*(sin(a*l_star*(t-t1_star))*cos(t)+a*l_star*cos(a*l_star*(t-t1_star))*sin(t))*(sin(i*l_star*(t-t1_star))*sin(t)+(i*i)*(l_star*l_star)*sin(i*l_star*(t-t1_star))*sin(t)-i*l_star*cos(i*l_star*(t-t1_star))*cos(t)*2.0))/(b*d*sin(t))-(cos(d*j*(f-f1))*sin(b*d*(f-f1))*1.0/tan(t)*(sin(a*l_star*(t-t1_star))*cos(t)+a*l_star*cos(a*l_star*(t-t1_star))*sin(t))*(sin(i*l_star*(t-t1_star))*cos(t)+i*l_star*cos(i*l_star*(t-t1_star))*sin(t))*2.0)/(b*d*sin(t))));
							Z=FN*(d*y*sin(d*y*(f-f1))*sin(l_star*x*(t-t1_star))-(sin(d*y*(f-f1))*cos(t)*(sin(l_star*x*(t-t1_star))*cos(t)+l_star*x*cos(l_star*x*(t-t1_star))*sin(t)))/(d*y)+(sin(d*y*(f-f1))*sin(t)*(sin(l_star*x*(t-t1_star))*sin(t)+(l_star*l_star)*(x*x)*sin(l_star*x*(t-t1_star))*sin(t)-l_star*x*cos(l_star*x*(t-t1_star))*cos(t)*2.0))/(d*y))/sin(t);
							coef(ii)=Z; ii++;
						}
					}
				}
			}
		}
	}
	
	return coef;
}