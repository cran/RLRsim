#include <R.h>
#include <Rmath.h>
#include <iostream>

using namespace std;

extern "C"{
	void RLRsim (int *p, //no. of fixed effects
			int *k, //min(no. of obs, no. of random effects)
			int *n, //no. of obs
			int *s, //no of values to simulate
			int *g, //length of lambda grid
			int *q, //no of tested fixed effects
			const double *mu, //scaled eigenvalues
			const double *lambda, //grid of vlaues for lambda
			const double *lambda0, //value of lambda under H0
			double *xi, //scaled eigenvalues
			const bool *REML,
			double *res)
{
GetRNGstate();

int is, ig, ik;

double **lambdamu = new double*[*g];
double **lambdamuP1 = new double*[*g];
double **fN = new double*[*g];
double **fD = new double*[*g];
for(ig=0; ig < *g; ++ig){
	lambdamu[ig] = new double[*k];
	lambdamuP1[ig] = new double[*k];
	fN[ig] = new double[*k];
	fD[ig] = new double[*k];
};
//double *lambdamuP1 = new double[*g][*k];
//double *fN = new double[*g][*k];
//double *fD = new double[*g][*k];

double *sumlog1plambdaxi = new double[*g];
double *Chi1 = new double[*k];

double  LR, N, D, ChiK, ChiSum;

int dfChiK, n0;



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

for(ig=0; ig < *g; ++ig){
	delete lambdamu[ig];
	delete lambdamuP1[ig];
	delete fN[ig];
	delete fD[ig];
};
delete [] lambdamu;
delete [] lambdamuP1;
delete [] fN;
delete [] fD;
delete [] sumlog1plambdaxi;
delete [] Chi1;

PutRNGstate();
}
}
