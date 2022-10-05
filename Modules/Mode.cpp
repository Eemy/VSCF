#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <stdbool.h>
#include "Mode.h"
#include "../UtilityFunc/ghQuad.h"
#include "../UtilityFunc/aux.h"

double pi =  3.1415926535897932384626433832795;
double au_to_wn = 219474.6313708;
double mass_au  = 1822.8884848961380701;

Mode::Mode(double _omega, int _nPoints, int _conv) {
    omega = _omega/au_to_wn;
    nPoints = _nPoints;
    nBasis = nPoints-1;
    //alpha = m*omega;

    waveFcn = new double[nBasis];
    oldWaveFcn = NULL;
    points = new double[nPoints];
    weights = new double[nPoints];
    hermiteEval = new double[nBasis*nPoints];
    norm = new double[nBasis];
    groundState = NULL;
    excitedState = NULL;

    //Set up some of the arrays
    gauher(points,weights,nPoints);
    for(int i=0 ; i<nBasis ; i++) {
      for(int j=0 ; j<nPoints ; j++) {
        hermiteEval[i*nPoints+j] = hermite(i,points[j]);
      }
      norm[i] = 1/(sqrt(pow(2.0,i)*factorial(i)))*pow(1/pi,0.25);
    }

    //Dipoles: Mode is excited
    excited = false;
    
    conv = _conv;
    //DIIS Set-up
    if(conv == 2) {
      maxDiisError = 0.0;
      diis_subspace = 6;
      Fsave.resize(diis_subspace);
      Esave.resize(diis_subspace); 
      density = new double[nBasis*nBasis];
    }
}

Mode::~Mode() {
  delete[] waveFcn;
  delete[] oldWaveFcn;
  delete[] weights;
  delete[] points;
  delete[] hermiteEval;
  delete[] norm;
  //delete[] groundState;
  delete[] excitedState;
  for(int i=0 ; i<diis_subspace ; i++) {
    delete[] Fsave[i];
    delete[] Esave[i];
  }
  if(conv == 2)
    delete[] density;
}

//===================================================================
void Mode::setGroundState() {
  groundState = new double[nBasis];
  for(int i=0 ; i<nBasis; i++)
    groundState[i] = waveFcn[i]; 
}

void Mode::setExcitedState() {
  excitedState = new double[nBasis];
  for(int i=0 ; i<nBasis; i++)
    excitedState[i] = waveFcn[i]; 
}

void Mode::updateWaveFcn(double* newWaveFcn) {
  if(oldWaveFcn != NULL)
    delete[] oldWaveFcn; //deallocate original memory block
  oldWaveFcn = waveFcn; //move pointer to another memory block
  waveFcn = newWaveFcn; 
  if(conv == 2)
    updateDensity();
}

void Mode::updateDensity() {
  for(int i=0 ; i<nBasis ; i++) {
    for(int j=0 ; j<nBasis ; j++) {
      density[i*nBasis+j] = waveFcn[i]*waveFcn[j];
    }
  }
}

double Mode::computeMaxDiff() {
  double diff = 0.0;
  for(int i=0 ; i<nBasis ; i++) {
    double temp = fabs(fabs(oldWaveFcn[i])-fabs(waveFcn[i]));
    if(temp > diff)
      diff = temp;
  }
  return diff;
}

//============================GETTERS/SETTERS================================
double* Mode::getWaveFcn() {return waveFcn;}
double* Mode::getGState() {return groundState;}
double* Mode::getEState() {return excitedState;}
double Mode::getOverlapEG() {
  double integral = 0.0;
  for(int i=0 ; i<nBasis ; i++)
    integral += groundState[i]*excitedState[i];
  return integral;
}
//double Mode::getAlpha() {return alpha;}
double Mode::getOmega() {return omega;}
//double Mode::getMass() {return m;}
int Mode::getNPoints() {return nPoints;}
double Mode::getWeight(int index) {return weights[index];}
double Mode::getHerm(int herm, int point) {return hermiteEval[herm*nPoints+point];}
double Mode::getNorm(int index) {return norm[index];}
double Mode::getDIISError() {return maxDiisError;}
double* Mode::getDensity() {return density;}

//double Mode::getPoint(int index) {return points[index]/(sqrt(alpha));}

void Mode::setExcited(bool status) {excited = status;}
//====================FOR EFFECTIVE POTENTIAL=======================
double Mode::getIntegralComponent(int point) {
  //Dipole: bra(excited) and ket(ground) different
  if(excited) {
    double groundIntegral = 0.0;
    double excitedIntegral = 0.0;
    for(int i=0 ; i<nBasis ; i++) {
      groundIntegral += hermiteEval[i*nPoints+point]*norm[i]*groundState[i];
      excitedIntegral += hermiteEval[i*nPoints+point]*norm[i]*excitedState[i];
    }
    return excitedIntegral*groundIntegral;
  }
  double integralComponent = 0.0;
  for(int i=0 ; i<nBasis ; i++) {
    integralComponent += hermiteEval[i*nPoints+point]*norm[i]*waveFcn[i];
  }

  return integralComponent*integralComponent;
}
//========================================DIIS Shenanigans=========================================
void Mode::diis(double *F, double *E, int iter) {
  printf("Iteration: %d\n",iter);  
  //Make copy of current Fock Matrix to save, the one passed in can be modified
  double *Fcopy = new double[nBasis*nBasis];
  std::copy(F,F+(nBasis*nBasis),Fcopy);
  setMaxElement(E);

  int index;
  if(iter>=0 && iter<diis_subspace) index = iter+1;
  if(iter>=diis_subspace) { 
    index = diis_subspace;
    delete[] Fsave[iter%diis_subspace];
    delete[] Esave[iter%diis_subspace];
  }
    Fsave[iter%diis_subspace] = Fcopy;
    Esave[iter%diis_subspace] = E; 

  //Extrapolate the Fock out of this thing -Justin
  if(iter>1) { //why not iter>0?
    //set up matrix [B11 B12 B13 ... -1.0]
    //              [B21 B22 B23 ... -1.0]
    //              [....................]
    //              [-1.0-1.0-1.0 ... 0.0]
    double *A = new double[(index+1)*(index+1)];
    for(int i=0 ; i<index+1 ; i++) A[i*(index+1)+index] = -1.0;
    for(int i=0 ; i<index+1 ; i++) A[index*(index+1)+i] = -1.0;
    A[index*(index+1)+index] = 0.0;

    for(int i=0 ; i<index ; i++) {
      for(int j=0 ; j<index ; j++) {
        A[i*(index+1)+j] = 0.0;
        for(int k=0 ; k<nBasis*nBasis ; k++) {
          A[i*(index+1)+j] += Esave[i][k]*Esave[j][k];//Bij
        }
      }
    }
    printmat(A,index+1,index+1);
     
    //solve for regression coeff
    double *B = new double[index+1];
    for(int i=0 ; i<index+1 ; i++) B[i] = 0.0;
    B[index] = -1.0;
    linsolver(A,B,index+1);

   //use coefficients for new Fock matrix
    for(int i=0 ; i<nBasis*nBasis ; i++) F[i] = 0.0;
    for(int i=0 ; i<index ; i++) {
      for(int j=0 ; j<nBasis*nBasis ; j++) {
        F[j] += B[i]*Fsave[i][j];
      }
    }
    delete[] A;
    delete[] B;
  }
}

void Mode::setMaxElement(double *array) {
  double max = 0.0;
  for(int i=0 ; i<nBasis*nBasis ; i++) {
    if (max<array[i]) max=array[i];
  }
  maxDiisError = max;
}
