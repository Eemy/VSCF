#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include "Mode.h"
#include "EigSolver.h"
#include "../UtilityFunc/aux.h"

EigSolver::EigSolver(int _nPoints, int _conv) {
  nPoints = _nPoints;
  nBasis = nPoints-1;
  T = new double[(nBasis+2)*(nBasis+2)];
  buildTmat();

  //conv: 1=roothaan 2=diis
  conv = _conv;
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
double EigSolver::solveMode(Mode* mode, std::vector<double> pot, int state, int iter) {
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
    printf("Fock Matrix\n"); 
    printmat(H,1,nBasis,nBasis,1.0);    
 
    //Compute Error Vector
    if(iter >= 0 && conv == 2)
      mode->saveErrorVec(H,iter);

    double* evals = new double[nBasis];
    diagonalize(H,evals,nBasis); //Hamiltonian becomes eigvec matrix
  
    //Update and Return Requested Information
    double evalNeeded = evals[state];
    mode->updateAllPsi_AllE(H,evals);

    
/*    //Compute Error Vector
    if(iter >= 0 && conv == 2)
      mode->saveErrorVec(iter);
*/
/*
    for(int i=0 ; i<nBasis ; i++) {
      printf("EVALS: %.8f\n",evals[i]*219474.6313708);
    }
*/
    //printf("EigVal:\n%.8f\n",evalNeeded*219474.6313708);

//    delete[] H;
//    delete[] evals;
    return evalNeeded;
}
//===============================================================
//=======================DIIS Shenanigans========================
void EigSolver::diis(std::vector<Mode*> dof) {
//void EigSolver::diis(Mode* mode) {
/*
  printf("Iteration: %d\n",iter);  
  setMaxElement(E);
*/
  for(int a=0 ; a<dof.size() ; a++) {
  int index = dof[a]->getNumErrorVecs();
  //int index = mode->getNumErrorVecs();

  //Extrapolate the Fock out of this thing -Justin
  if(index>=2) { //why not iter>0?
    //set up matrix [B11 B12 B13 ... -1.0]
    //              [B21 B22 B23 ... -1.0]
    //              [....................]
    //              [-1.0-1.0-1.0 ... 0.0]
    double *A = new double[(index+1)*(index+1)];
    for(int i=0 ; i<(index+1)*(index+1) ; i++) A[i] = 0.0;
    for(int i=0 ; i<index+1 ; i++) A[i*(index+1)+index] = -1.0;
    for(int i=0 ; i<index+1 ; i++) A[index*(index+1)+i] = -1.0;
    A[index*(index+1)+index] = 0.0;

    //for(int a=0 ; a<dof.size() ; a++) {
      for(int i=0 ; i<index ; i++) {
        for(int j=0 ; j<index ; j++) {
          A[i*(index+1)+j] += dof[a]->dotErrorVecs(i,j);
//          A[i*(index+1)+j] += mode->dotErrorVecs(i,j);
        }
      }
    //}
    printmat(A,1,index+1,index+1,1.0);
     
    //solve for regression coeff
    double *B = new double[index+1];
    for(int i=0 ; i<index+1 ; i++) B[i] = 0.0;
    B[index] = -1.0;
    linsolver(A,B,index+1);

    printmat(B,1,1,index+1,1.0);

    //use coefficients for new density matrix
    //for(int a=0 ; a<dof.size() ; a++) {
      dof[a]->extrapolateDensity(B);
    //}

//    mode->extrapolateDensity(B);
    //for(int i=0 ; i<nBasis*nBasis ; i++) F[i] = 0.0;
    //for(int i=0 ; i<index ; i++) {
    //  for(int j=0 ; j<nBasis*nBasis ; j++) {
    //    F[j] += B[i]*Fsave[i][j];
    //  }
    //}
    delete[] A;
    delete[] B;
  }
  }
}
