#include <stdio.h>
#include <cmath>
#include "ghQuad.h"

void gauher(double *x, double *w, int n)
{
	const double EPS=1.0e-14;
  const double PIM4=0.7511255444649425;
	const int MAXIT=100;
	int i,its,j,m;
  double p1,p2,p3,pp,z,z1;

	m=(n+1)/2;
	for (i=0;i<m;i++) {
		if (i == 0) {
			z=sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
		} else if (i == 1) {
			z -= 1.14*pow((double)(n),0.426)/z;
		} else if (i == 2) {
			z=1.86*z-0.86*x[0];
		} else if (i == 3) {
			z=1.91*z-0.91*x[1];
		} else {
			z=2.0*z-x[i-2];
		}
		for (its=0;its<MAXIT;its++) {
			p1=PIM4;
			p2=0.0;
			for (j=0;j<n;j++) {
				p3=p2;
				p2=p1;
				p1=z*sqrt(2.0/(j+1))*p2-sqrt((double)(j)/(j+1))*p3;
			}
			pp=sqrt((double)(2*n))*p2;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its >= MAXIT) printf("too many iterations in gauher");
		x[i]=z;
		x[n-1-i] = -z;
		w[i]=2.0/(pp*pp);
		w[n-1-i]=w[i];
	}
}

double factorial(int n)
{
  double fac = 1.0;
  if(n>1)
  {
    for(int i=2 ; i<=n ; i++)
      fac *= (double)i;
  }  
  return fac;
}

//***********************************************
double hermite(int n, double x)
{
  double herm;

  if(n==0)
    herm = 1.0;
  else if(n==1)
    herm = 2.0*x;
  else if(n==2)
    herm = 4.0*x*x-2.0;
  else 
    herm = 2.0*x*hermite(n-1,x)-2.0*((double)(n-1))*hermite(n-2,x);
 
  return herm;
 
}  
