#include "spectdg.h"

void SpectDG::spatialMat()
//system ODE matrices are built by integrating weighted parts of the governing equations computed in all points
{
	int ii=0, i, j, a, b, x, y;
	double V1, V2, V3, V4, V5, V6, FN, Z, tau_phi;
	ncoef = ncoef1;
	vec coef1 = quad2d(&SpectDG::coefFun1, t1, t2, f1, f2);
	ncoef = ncoef2;
	vec coef2 = quad2d(&SpectDG::coefFun2, t1_star, t2_star, f1, f2);
	
	BN=zeros(M*N+M1*N1);

	//large-scale
	for(j=1;j<=2*N-1;j+=2){
		for(i=1;i<=M;i++){
			h1(i - 1 + (j - 1)*M / 2) =coef1[ii] ; ii++;
			h2(i - 1 + (j - 1)*M / 2) = coef1[ii] ; ii++;
			h3(i - 1 + (j - 1)*M / 2) =coef1[ii] ; ii++;
			BW(i-1+(j-1)*M/2)=coef1[ii] ; ii++;
			for(b=1;b<=2*N-1;b+=2){
				for(a=1;a<=M;a++){
					A(i - 1 + (j - 1)*M / 2, a - 1 + (b - 1)*M / 2) = coef1[ii] ; ii++;
					H(i - 1 + (j - 1)*M / 2, a - 1 + (b - 1)*M / 2) =  coef1[ii] ; ii++;
					JJ(i - 1 + (j - 1)*M / 2, a - 1 + (b - 1)*M / 2) = coef1[ii] ; ii++;
					C(i - 1 + (j - 1)*M / 2, a - 1 + (b - 1)*M / 2) = coef1[ii] ; ii++;
					J(i - 1 + (j - 1)*M / 2, a - 1 + (b - 1)*M / 2) = coef1[ii] ; ii++;
					D(i - 1 + (j - 1)*M / 2, a - 1 + (b - 1)*M / 2) = coef1[ii] ; ii++;
				}
			}
		}	
	}

	//convective terms
	
	for(y=1;y<=2*N-1;y+=2){
		for(x=1;x<=M;x++){
			for(b=1;b<=2*N-1;b+=2){
				for(a=1;a<=M;a++){
					for(j=1;j<=2*N-1;j+=2){
						for(i=1;i<=M;i++){
							BNC((a-1+M*(b-1)/2),(i-1+M*(j-1)/2),(x-1+(y-1)*M/2))=coef1(ii); ii++;
						}
					}
				}
			}
		}
	}

	ii = 0 ;

	//small-scale

	for(j=1;j<=2*N1-1;j+=2){
		for(i=1;i<=M1;i++){
			h1(M*N+i - 1 + (j - 1)*M1 / 2) = coef2[ii] ; ii++;
			h2(M*N+i - 1 + (j - 1)*M1 / 2) = coef2[ii] ; ii++;
			h3(M*N+i - 1 + (j - 1)*M1 / 2) = coef2[ii] ; ii++;
			BW(M*N+i-1+(j-1)*M1/2)= coef2[ii] ; ii++;
			for(b=1;b<=2*N1-1;b+=2){
				for(a=1;a<=M1;a++){
					A(M*N+i - 1 + (j - 1)*M1 / 2, M*N+a - 1 + (b - 1)*M1 / 2) =  coef2[ii] ; ii++;
					H(M*N+i - 1 + (j - 1)*M1 / 2, M*N+a - 1 + (b - 1)*M1 / 2) = coef2[ii] ; ii++;
					JJ(M*N+i - 1 + (j - 1)*M1 / 2, M*N+a - 1 + (b - 1)*M1 / 2) = 0.0; ii++;//coef2[ii] ; ii++;
					C(M*N+i - 1 + (j - 1)*M1 / 2, M*N+a - 1 + (b - 1)*M1 / 2) =  coef2[ii] ; ii++;
					J(M*N+i - 1 + (j - 1)*M1 / 2, M*N+a - 1 + (b - 1)*M1 / 2) = 0.0; ii++;//coef2[ii] ; ii++;
					D(M*N+i - 1 + (j - 1)*M1 / 2, M*N+a - 1 + (b - 1)*M1 / 2) = coef2[ii] ; ii++;
				}
			}
		}
	}
	
	//convective terms

	for(y=1;y<=2*N1-1;y+=2){
		for(x=1;x<=M1;x++){
			for(b=1;b<=2*N1-1;b+=2){
				for(a=1;a<=M1;a++){
					for(j=1;j<=2*N1-1;j+=2){
						for(i=1;i<=M1;i++){
							BNC((M*N+a-1+M1*(b-1)/2),(M*N+i-1+M1*(j-1)/2),(M*N+x-1+(y-1)*M1/2))=coef2(ii); ii++;
						}
					}
				}
			}
		}
	}
}