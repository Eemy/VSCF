#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "Mode.h"
#include "../UtilityFunc/ghQuad.h"

double pi =  3.1415926535897932384626433832795;
double au_to_wn = 219474.6313708;
double mass_au  = 1822.8884848961380701;

Mode::Mode(double _omega, double _m, int _nPoints) {
    omega = _omega/au_to_wn;
    m = _m*mass_au;
    nPoints = _nPoints;
    nBasis = nPoints-1;
    alpha = m*omega;

    waveFcn = new double[nBasis];
    oldWaveFcn = NULL;
    points = new double[nPoints];
    weights = new double[nPoints];
    hermiteEval = new double[nBasis*nPoints];
    norm = new double[nBasis];
    groundState = NULL;

    //Set up some of the arrays
    gauher(points,weights,nPoints);
    for(int i=0 ; i<nBasis ; i++) {
      for(int j=0 ; j<nPoints ; j++) {
        hermiteEval[i*nPoints+j] = hermite(i,points[j]);
      }
      norm[i] = 1/(sqrt(pow(2.0,i)*factorial(i)))*pow(1/pi,0.25);
    }
}

Mode::~Mode() {
  delete[] waveFcn;
  delete[] oldWaveFcn;
  delete[] weights;
  delete[] points;
  delete[] hermiteEval;
  delete[] norm;
  delete[] groundState;
}

//===================================================================
void Mode::setGroundState() {
  groundState = waveFcn; 
}

void Mode::updateWaveFcn(double* newWaveFcn) {
  delete[] oldWaveFcn; //deallocate original memory block
  oldWaveFcn = waveFcn; //move pointer to another memory block
  waveFcn = newWaveFcn; 
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

//============================GETTERS================================
double* Mode::getWaveFcn() {return waveFcn;}
double Mode::getAlpha() {return alpha;}
double Mode::getOmega() {return omega;}
double Mode::getMass() {return m;}
int Mode::getNPoints() {return nPoints;}
double Mode::getWeight(int index) {return weights[index];}
double Mode::getHerm(int herm, int point) {return hermiteEval[herm*nPoints+point];}
double Mode::getNorm(int index) {return norm[index];}
double Mode::getPoint(int index) {return points[index];}

//====================FOR EFFECTIVE POTENTIAL=======================
double Mode::getIntegralComponent(int point) {
  double integralComponent = 0.0;
  for(int i=0 ; i<nBasis ; i++) {
    integralComponent += hermiteEval[i*nPoints+point]*norm[i]*waveFcn[i];
  }
//  integralComponent *= weights[point];//deleted sqrt(alpha);
  return integralComponent;
}
