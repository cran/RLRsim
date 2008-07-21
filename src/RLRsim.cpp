#include <R.h>
#include <Rmath.h>
#include <iostream>

using namespace std;

extern "C"{
void RLRsim (const int *p,
             const int *k,
													const int *n,
													const int *s,
													const int *g,
													const int *q,
													const double *mu,
													const double *lambda,
													const double *lambda0,
													double *xi,
													const bool *REML,
													double *res)
{
GetRNGstate();


double lambdamu[*g][*k];
double lambdamuP1[*g][*k];
double fN[*g][*k];
double fD[*g][*k];
double sumlog1plambdaxi[*g];
double Chi1[*k];

double  LR, N, D, Chiq, ChiK, ChiSum;

int dfChiK, n0;

unsigned int is, ig, ik;

dfChiK = imax2(0, *n-*p-*k);

if(*REML == true){
		n0 = *n - *p;
		for(ik=0; ik < *k; ++ik)
		{
		  xi[ik] = mu[ik];
		}  
}	else {
	 n0 = *n;
}	 



/*precompute stuff that stays constant over simulations*/  
for(ig=0; ig< *g; ++ig)
{
	sumlog1plambdaxi[ig] = 0;
	
	for(ik=0; ik < *k; ++ik)
	{
		 lambdamu[ig][ik] = lambda[ig]*mu[ik];
	 	lambdamuP1[ig][ik] = lambdamu[ig][ik]+1.0;
		 
			fN[ig][ik] = ((lambda[ig] - (*lambda0))*mu[ik])/lambdamuP1[ig][ik]; 
			fD[ig][ik] = (1 + ((*lambda0)*mu[ik]) )/lambdamuP1[ig][ik];

		 sumlog1plambdaxi[ig] += log1p(lambda[ig]*xi[ik]);		
			 
	}
}


for(is=0; is< *s; ++is)
{

 /*make random variates, set LR 0*/
 LR =  0;
 ChiSum = 0;
 ChiK = rchisq(dfChiK);
 for(ik=0; ik < *k; ++ik)
	{
			Chi1[ik]=rchisq(1);
		 if(REML == false){
					 	ChiSum += Chi1[ik];
			} 	
 } 
 
 
	for(ig=1; ig < *g; ++ig)   /*loop over lambda-grid*/
	{
			N = D = 0;
	
			for(ik=0; ik < *k; ++ik)  /*loop over mu,xi*/
			{
						N += fN[ig][ik] * Chi1[ik];
						D += fD[ig][ik] * Chi1[ik];
  	}
  	
  	D = D + ChiK;

  	LR = n0 * log1p(N/D) -  sumlog1plambdaxi[ig];
			
			if(LR >= res[is]){   /*save if LR is bigger than previous LR*/
			 	res[is]= LR;
			}
			else break;         /*else keep larger value and go to next iteration */
			
 }/*end for *g*/ 
 
 /* add additional term for LR*/
	if(REML==false)
 {
  	 res[is] += (*n) * (log1p(rchisq(*q)/(ChiSum+ChiK)));
 }
}/*end for s*/ 

PutRNGstate();
}
}
