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

    //Hold modals (VMP2)
    energies = NULL; 
    waveAll = NULL;
    waveAll_prev = NULL;
    tempAll = NULL;
    //Hold VSCF states
    vscfPsi = new double[nBasis*2];

    //Control which states to use
    vscfStates = false; //default to modals
    harm = false;
    bra = 0;
    ket = 0;
    
    //Set up gauher and norm arrays 
    points = new double[nPoints];
    weights = new double[nPoints];
    hermiteEval = new double[nBasis*nPoints];
    norm = new double[nBasis];

    gauher(points,weights,nPoints);
    for(int i=0 ; i<nBasis ; i++) {
      for(int j=0 ; j<nPoints ; j++) {
        hermiteEval[i*nPoints+j] = hermite(i,points[j]);
      }
      norm[i] = 1/(sqrt(pow(2.0,i)*factorial(i)))*pow(1/pi,0.25);
    }
    density = new double[nBasis*nBasis];

    //DIIS Set-up
    conv = _conv;
    if(conv == 2) {
      diis_subspace = 5;
//      Fsave.resize(diis_subspace);
//      Dsave.resize(diis_subspace);
//      Esave.resize(diis_subspace); 
    }
}

Mode::~Mode() {
  delete[] weights;
  delete[] points;
  delete[] hermiteEval;
  delete[] norm;
  delete[] vscfPsi;
  delete[] waveAll;
  delete[] waveAll_prev;
  delete[] energies;
  delete[] density;
  if(conv == 2) {
    for(int i=0 ; i<Dsave.size() ; i++) {
//      delete[] Fsave[i];
      delete[] Dsave[i];
      delete[] Esave[i];
    }
  }
}

//===================================================================
void Mode::setGroundState() {
  for(int i=0 ; i<nBasis; i++)
    vscfPsi[i] = waveAll[i]; 
}

void Mode::setExcitedState() {
  for(int i=0 ; i<nBasis; i++)
    vscfPsi[nBasis+i] = waveAll[nBasis+i]; 
}

void Mode::updateAllPsi_AllE(double* newPsi, double* newE) {
  if(waveAll_prev != NULL)
    delete[] waveAll_prev; 
  if(waveAll != NULL)
    waveAll_prev = waveAll; //pass pointer of original memory block
  if(energies != NULL)
    delete[] energies; //deallocate original memory block
  waveAll = newPsi; //move pointer to another memory block
  energies = newE; 

//  updateDensity(); 
}

void Mode::updateDensity() {
  for(int i=0 ; i<nBasis ; i++) {
    for(int j=0 ; j<nBasis ; j++) {
      if(vscfStates)
        density[i*nBasis+j] = vscfPsi[bra*nBasis+i]*vscfPsi[ket*nBasis+j];
      else       
        density[i*nBasis+j] = waveAll[bra*nBasis+i]*waveAll[ket*nBasis+j];
    }
  }
}

//void Mode::extrapolateDensity(double *coeff) {
void Mode::extrapolateDensity(double *coeff,int iter) {
  if(conv==2) { //coeff are the LC coefficients 
    for(int i=0 ; i<nBasis*nBasis ; i++) density[i] = 0.0;
    for(int i=0 ; i<Esave.size() ; i++) {
      for(int j=0 ; j<nBasis*nBasis ; j++) {
        density[j] += coeff[i]*Dsave[i][j];
      }
    }
    /*//Replace saved density
    double *Dcopy = new double[nBasis*nBasis];
    std::copy(density,density+(nBasis*nBasis),Dcopy);
    delete [] Dsave[iter%diis_subspace];
    Dsave[iter%diis_subspace] = Dcopy;*/
  } else if(conv ==3) { //coeff is the Fock Matrix
    double learningFactor = 0.1;
    for(int i=0 ; i<nBasis*nBasis ; i++) {
      density[i] = density[i]+learningFactor*coeff[i];
    }
  }
}

void Mode::saveCurrentDensity(int iter) {
/*
  double *Dcopy = new double[nBasis*nBasis];
  std::copy(density,density+(nBasis*nBasis),Dcopy);

  if(iter==-1) {
    firstDensity = Dcopy; 
  } else {
*/
    double *Dcopy = new double[nBasis*nBasis];
    std::copy(density,density+(nBasis*nBasis),Dcopy);

    int index = iter;
    //int index = iter+1;

    if(index >=diis_subspace) { 
      delete[] Dsave[(index)%diis_subspace];
      Dsave[(index)%diis_subspace] = Dcopy;
    } else {
      Dsave.push_back(Dcopy);
    }
//  }
}

void Mode::saveErrorVec(double *F, int iter) {
//void Mode::saveErrorVec(int iter) {
/*  //Make copy of current Fock Matrix to save
  double *Fcopy = new double[nBasis*nBasis];
  std::copy(F,F+(nBasis*nBasis),Fcopy);*/

  //Compute Error Vector
  double *fd = new double[nBasis*nBasis];
  double *df = new double[nBasis*nBasis]; 
  double *error = new double[nBasis*nBasis];
  ABmult(fd,F,density,nBasis,nBasis,nBasis,nBasis,nBasis,nBasis,1);
  ABmult(df,density,F,nBasis,nBasis,nBasis,nBasis,nBasis,nBasis,1);
  for(int i=0 ; i<nBasis*nBasis ; i++) error[i] = fd[i]-df[i]; 
/*
  double *error = new double[nBasis*nBasis];
*/
/*
  if(iter==0) {
    for(int i=0 ; i<nBasis*nBasis ; i++) 
      error[i] = density[i]-firstDensity[i];
  } else {
*/
/*    for(int i=0 ; i<nBasis*nBasis ; i++) 
      error[i] = density[i]-Dsave[iter%diis_subspace][i];
*/
//      error[i] = density[i]-Dsave[(iter-1)%diis_subspace][i];
//  }

  printf("Error Matrix\n");
  printmat(error,1,nBasis,nBasis,1.0);
  setMaxElement(error);  
  if(iter >=diis_subspace) { 
    delete[] Esave[iter%diis_subspace];
    Esave[iter%diis_subspace] = error; 
  } else {
    Esave.push_back(error);
  }
  saveCurrentDensity(iter);
////debug
  printf("EigVecs\n");
  printmat(waveAll,1,nBasis,nBasis,1.0);
  printf("Density Matrix\n");
  printmat(density,1,nBasis,nBasis,1.0);
//  printf("Fock Matrix\n");
//  printmat(F,1,nBasis,nBasis,1.0);
////debug
  //Save Fock Matrix and Error Matrix
//  if(iter >=diis_subspace) { 
//    delete[] Fsave[iter%diis_subspace];
//    Fsave[iter%diis_subspace] = Fcopy;
//  }
   
//  delete[] fd;
//  delete[] df;
}

double Mode::dotErrorVecs(int i, int j) {
  double sum = 0.0;
  for(int k=0 ; k<nBasis*nBasis ; k++) {
    sum += Esave[i][k]*Esave[j][k];
  }
  return sum;
}

void Mode::resetSubspace() {
  for(int i=0 ; i<Dsave.size() ; i++) {
    delete[] Dsave[i];
  }
  for(int i=0 ; i<Esave.size() ; i++) {
    delete[] Esave[i];
  }
  Dsave.clear();
  Esave.clear();
  delete[] firstDensity;
}

void Mode::setMaxElement(double *array) {
  double max = 0.0;
  for(int i=0 ; i<nBasis*nBasis ; i++) {
    if (max<array[i]) max=array[i];
  }
  maxDiisError = max;
}

double Mode::computeMaxDiff() {
  double diff = 0.0;
  for(int i=0 ; i<nBasis ; i++) {
    double temp = fabs(fabs(waveAll_prev[bra*nBasis+i])-fabs(waveAll[bra*nBasis+i]));
    if(temp > diff)
      diff = temp;
  }
  return diff;
}

//============================GETTERS/SETTERS================================
double Mode::getOverlapEG() {
  double integral = 0.0;
  if(vscfStates) {
    for(int i=0 ; i<nBasis ; i++)
      integral += vscfPsi[i]*vscfPsi[nBasis+i];
  } else {
    for(int i=0 ; i<nBasis ; i++)
      integral += waveAll[bra*nBasis+i]*waveAll[ket*nBasis+i];
  }
  return integral;
}

double* Mode::getModal(int state) {
  double* modal = new double[nBasis];
  for(int i=0 ; i<nBasis ; i++) 
    modal[i] = waveAll[nBasis*state+i];
  return modal;
}
double Mode::getEModal(int state) {return energies[state];}
double Mode::getEModal() {return energies[ket];}
int Mode::getNumErrorVecs() {return Esave.size();}
double Mode::getDIISError() {return maxDiisError;}

//double Mode::getAlpha() {return alpha;}
//double Mode::getMass() {return m;}
//double Mode::getPoint(int index) {return points[index]/(sqrt(alpha));}
double Mode::getOmega() {return omega;}
int Mode::getNPoints() {return nPoints;}
double Mode::getWeight(int index) {return weights[index];}
double Mode::getHerm(int herm, int point) {return hermiteEval[herm*nPoints+point];}
double Mode::getNorm(int index) {return norm[index];}
double* Mode::getDensity() {return density;}
void Mode::setBra(int _bra) {
  if(vscfStates && _bra>1) {
    printf("Error: No existing VSCF states above 1 quanta.\n"); 
    exit(0);
  }
  bra = _bra; 
  updateDensity();
}
void Mode::setKet(int _ket) {
  if(vscfStates && _ket>1) {
    printf("Error: No existing VSCF states above 1 quanta.\n"); 
    exit(0);
  }
  ket = _ket; 
  updateDensity();
}
int Mode::getBra() {return bra;}
int Mode::getKet() {return ket;}
void Mode::useVSCFStates(bool use) {
  vscfStates = use;
  if(vscfStates && (bra>1 || ket>1)) {
    printf("Error: No existing VSCF states above 1 quanta. Change bra and ket and try again.\n"); 
    vscfStates = false;
    exit(0);
  }
  updateDensity(); 
}

void Mode::setHarmonic() {
  double* harmPsi = new double[nBasis*nBasis];
  for(int i=0 ; i<nBasis ; i++) {
    for(int j=0 ; j<nBasis ; j++) {
      if(i==j)
        harmPsi[i*nBasis+j] = 1;
      else
        harmPsi[i*nBasis+j] = 0;
    }
  }
  if(waveAll != NULL)
    tempAll = waveAll; //pass pointer of anharmonic states, if present 
  waveAll = harmPsi;
  harm = true;
}

void Mode::setAnharmonic() {
  if(tempAll != NULL && harm) {
    delete[] waveAll;
    waveAll = tempAll;
    tempAll = NULL; 
    harm = false;
  } else {
    printf("setHarmonic() was never called or there are no available anharmonic states.\n");
  }
}
//====================FOR EFFECTIVE POTENTIAL======================
double Mode::getIntegralComponent(int point) {
  double integralComponent = 0.0;
/*
  if(vscfStates) {
    if(bra > 1 || ket > 1) {
      printf("Error: No existing VSCF states above 1 quanta of excitation.\n"); 
      exit(0);
    }
    double braIntegral = 0.0;
    double ketIntegral = 0.0;
    for(int i=0 ; i<nBasis ; i++) {
      braIntegral += hermiteEval[i*nPoints+point]*norm[i]*vscfPsi[bra*nBasis+i];
      ketIntegral += hermiteEval[i*nPoints+point]*norm[i]*vscfPsi[ket*nBasis+i];
    }
    return braIntegral*ketIntegral*weights[point];
  }
  //Else use modal states
  double braIntegral = 0.0;
  double ketIntegral = 0.0;
  for(int i=0 ; i<nBasis ; i++) {
    braIntegral += hermiteEval[i*nPoints+point]*norm[i]*waveAll[bra*nBasis+i];
    ketIntegral += hermiteEval[i*nPoints+point]*norm[i]*waveAll[ket*nBasis+i];
  }
  return braIntegral*ketIntegral*weights[point];
*/

  for(int i=0 ; i<nBasis ; i++) {
    for(int j=0 ; j<nBasis ; j++) {
      integralComponent += density[i*nBasis+j]*hermiteEval[i*nPoints+point]*norm[i]*hermiteEval[j*nPoints+point]*norm[j];
    }
  }
  return integralComponent*weights[point];

}

