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

    //DIIS Set-up
    conv = _conv;
    if(conv == 2) {
      maxDiisError = 0.0;
      diis_subspace = 5;
      Fsave.resize(diis_subspace);
      Esave.resize(diis_subspace); 
      density = new double[nBasis*nBasis];
    }
}

Mode::~Mode() {
//  delete[] waveFcn;
//  delete[] oldWaveFcn;
  delete[] weights;
  delete[] points;
  delete[] hermiteEval;
  delete[] norm;
//  delete[] groundState;
//  delete[] excitedState;
  delete[] vscfPsi;
  delete[] waveAll;
  delete[] waveAll_prev;
  delete[] energies;
  for(int i=0 ; i<diis_subspace ; i++) {
    delete[] Fsave[i];
    delete[] Esave[i];
  }
  if(conv == 2)
    delete[] density;
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

  if(conv == 2) 
    updateDensity(); 
}

void Mode::updateDensity() {
  for(int i=0 ; i<nBasis ; i++) {
    for(int j=0 ; j<nBasis ; j++) {
      density[i*nBasis+j] = waveAll[bra*nBasis+i]*waveAll[bra*nBasis+j];
    }
  }
}

void Mode::saveErrorVec(double *F, int iter) {
  //Make copy of current Fock Matrix to save
  double *Fcopy = new double[nBasis*nBasis];
  std::copy(F,F+(nBasis*nBasis),Fcopy);

  //Compute Error Vector
  double *fd = new double[nBasis*nBasis];
  double *df = new double[nBasis*nBasis]; 
  double *error = new double[nBasis*nBasis];
  ABmult(fd,F,density,nBasis,nBasis,nBasis,nBasis,nBasis,nBasis,1);
  ABmult(df,density,F,nBasis,nBasis,nBasis,nBasis,nBasis,nBasis,1);
  for(int i=0 ; i<nBasis*nBasis ; i++) error[i] = fd[i]-df[i]; 

  //Save Fock Matrix and Error Matrix
  if(iter>=diis_subspace) { 
    delete[] Fsave[iter%diis_subspace];
    delete[] Esave[iter%diis_subspace];
  }
  Fsave[iter%diis_subspace] = Fcopy;
  Esave[iter%diis_subspace] = error; 
   
  delete[] fd;
  delete[] df;
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

//double Mode::getAlpha() {return alpha;}
//double Mode::getMass() {return m;}
//double Mode::getPoint(int index) {return points[index]/(sqrt(alpha));}
double Mode::getOmega() {return omega;}
int Mode::getNPoints() {return nPoints;}
double Mode::getWeight(int index) {return weights[index];}
double Mode::getHerm(int herm, int point) {return hermiteEval[herm*nPoints+point];}
double Mode::getNorm(int index) {return norm[index];}
double Mode::getDIISError() {return maxDiisError;}
double* Mode::getDensity() {return density;}
void Mode::setBra(int _bra) {bra = _bra;}
void Mode::setKet(int _ket) {ket = _ket;}
int Mode::getBra() {return bra;}
int Mode::getKet() {return ket;}
void Mode::useVSCFStates(bool use) {vscfStates = use;}
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
//====================FOR EFFECTIVE POTENTIAL=======================
double Mode::getIntegralComponent(int point) {
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
    return braIntegral*ketIntegral;
  }

  //Else use modal states
  double braIntegral = 0.0;
  double ketIntegral = 0.0;
  for(int i=0 ; i<nBasis ; i++) {
    braIntegral += hermiteEval[i*nPoints+point]*norm[i]*waveAll[bra*nBasis+i];
    ketIntegral += hermiteEval[i*nPoints+point]*norm[i]*waveAll[ket*nBasis+i];
  }
  return braIntegral*ketIntegral;
}

