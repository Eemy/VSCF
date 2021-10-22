#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "Mode.h"
#include "EigSolver.h"
#include "../UtilityFunc/aux.h"

void buildTmat();

EigSolver::EigSolver(int _nPoints) {
  nPoints = _nPoints;
  nBasis = nPoints-1;
  T = new double[(nBasis+2)*(nBasis+2)];
  buildTmat();
}

EigSolver::~EigSolver(){
  delete[] T;
}

//==============================================================
void EigSolver::buildTmat() {
  double* a = new double[(nBasis+2)*(nBasis+2)]();
  double* adag = new double[(nBasis+2)*(nBasis+2)]();
  for(int j=0 ; j<(nBasis+2) ; j++) {
    for(int k=0 ; k<(nBasis+2); k++) {
      if(j==(k-1))
        a[j*(nBasis+2)+k] = sqrt((double) k);
      if(j==(k+1))
        adag[j*(nBasis+2)+k] = sqrt((double) j);
    }
  }

  double* p = new double[(nBasis+2)*(nBasis+2)]();
  for(int i=0 ; i<(nBasis+2)*(nBasis+2) ; i++) {
    p[i] = sqrt(0.5)*(adag[i]-a[i]);//don't forget to add -0.25 prefactor with omega
  }
  ABmult(T,p,p,(nBasis+2),(nBasis+2),(nBasis+2),(nBasis+2),(nBasis+2),(nBasis+2),1);

  delete[] a;
  delete[] adag;
  delete[] p;
}

//===============================================================
double EigSolver::solveMode(Mode* mode, double* pot, int state) {
    double* H = new double[nBasis*nBasis];
    for(int i=0 ; i<nBasis ; i++) {
      for(int j=0 ; j<nBasis ; j++) {
        double VMat = 0.0;
        for(int k=0 ; k<nPoints ; k++) {
          VMat += mode->getHerm(i,k)*mode->getHerm(j,k)*mode->getNorm(i)*mode->getNorm(j)
                          *mode->getWeight(k)*pot[k];
        }
        H[i*nBasis+j] = mode->getOmega()*-0.5*T[i*(nBasis+2)+j]+VMat;
      }
    }
/*
    for(int i=0 ; i<nBasis ; i++) {
      for(int j=0 ; j<nBasis ; j++) {
        printf("H: %.8f\n",H[i*nBasis+j]);
      }
      printf("\n");
    }    
*/
    double* evals = new double[nBasis];
    diagonalize(H,evals,nBasis); //Hamiltonian becomes eigvec matrix
  
    //Update and Return Requested Information
    double* newWaveFcn = new double[nBasis];
    for(int i=0 ; i<nBasis ; i++) {
//      printf("%.8f\n",H[nBasis*state+i]);
      newWaveFcn[i] = H[nBasis*state+i];
    }
    mode->updateWaveFcn(newWaveFcn);
      double evalNeeded = evals[state];
/*
    for(int i=0 ; i<nBasis ; i++) {
      printf("EVALS: %.8f\n",evals[i]*219474.6313708);
    }
    for(int i=0 ; i<nBasis ; i++) {
      for(int j=0 ; j<nBasis ; j++) {
        printf("EVECS: %.8f\n",evecs[i*nBasis+j]);
      }
      printf("\n");
    }
*/
    //printf("EigVal:\n%.8f\n",evalNeeded*219474.6313708);

    delete[] H;
    delete[] evals;
//    delete[] evecs;
    return evalNeeded;
}
//===============================================================
