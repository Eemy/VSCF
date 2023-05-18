#include "Mp2Corr.h"
#include "Mode.h"
#include "Potential.h"
#include "../UtilityFunc/aux.h"
#include <vector>
#include <utility>
#include <iostream>
#include <cmath>
#include <stdbool.h>
using std::vector;
using std::pair;

Mp2Corr::Mp2Corr(vector<Mode*>& _dof, vector<Potential*>& _pot, vector<vector<vector<int>>>& _potIterators) { 
  dof = _dof;
  pot = _pot;
  tuples = _potIterators;

  maxState = 4;
  minState = 1;
  numStates = maxState-minState+1;
  nModes = dof.size();
  numPsi = 0;

  maxDim = 0;
  for(int i=0 ; i<pot.size() ; i++) {
    if(pot[i]->dim > maxDim)
      maxDim = pot[i]->dim;
  }
}

Mp2Corr::~Mp2Corr() {

}

void Mp2Corr::calculateIntegrals(vector<pair<vector<int>,vector<int>>>psi_m, vector<double> excitedEnergies) {
  //Set up
  numPsi = psi_m.size(); //size of degenerate subspace 

  singles.resize(numPsi*nModes*numStates);
  singlesDenom.resize(numPsi*nModes*numStates);
  for(int i=0 ; i<maxDim ; i++) {
    int excitationLevel = i+2;
    int size = numPsi;
    for(int j=0 ; j<excitationLevel ; j++) { //num unique tuples 
      size *= nModes-j; 
      size /= j+1;
    }
    size *= pow(numStates,excitationLevel);
    vector<double> integral(size);
    vector<double> denominator(size);
    integrals.push_back(integral);
    denominators.push_back(denominator);
  } 
  
  //Start calculation
  for(int i=0 ; i<numPsi ; i++) {
    for(int i2=0 ; i2<psi_m[i].first.size() ; i2++) { 
      int mode = psi_m[i].first[i2];
      int state = psi_m[i].second[i2]; 
      dof[mode]->setBra(state);
    }
    
    for(int j=0 ; j<tuples.size() ; j++) {
    for(int j2=0 ; j2<tuples[j].size() ; j2++) {

      //Single Excitations 
      for(int k=0 ; k<nModes ; k++) {
        vector<int> diff1 = psi_m[i].first;
        //diff1.push_back(i); //change for psi_m
        diff1.push_back(k);
        for(int l=minState ; l<=maxState ; l++) {
          dof[k]->setKet(l);
          if(integralIsNonZero(diff1,tuples[j][j2]) && !brillouin()) { 
            double integralVal = pot[j]->integrateTuple(j2,false);
            int index = i*numPsi*numStates+k*numStates+(l-minState);//change for psi_m
            //printf("Integral: %i %i State: %i, Tuple Num: %i Val: %.12f\n",i,k,l,j2,integralVal);   
            singles[index] += integralVal;
            if(singlesDenom[index] == 0.0)
              singlesDenom[index] = excitedEnergies[i+1]-(excitedEnergies[0]+(dof[k]->getEModal()-dof[k]->getEModal(0)));
          }

          //Double+ Excitations
          fillCorrectionMatrices(2, i, j, j2, k+1, diff1, excitedEnergies);
          dof[k]->setKet(0);
        }//l
      }//k
    }//j2
    }//j
    
    for(int i2=0 ; i2<psi_m[i].first.size() ; i2++) { 
      int mode = psi_m[i].first[i2];
      dof[mode]->setBra(0);
    }
  }//i

  printf("Singles\n");
  double *singlesVector = &singles[0];
  printmat(singlesVector,numPsi,nModes,numStates,1.0);
  printf("Doubles\n");
  double *doublesVector2 = &integrals[0][0];
  int nPairs = (nModes)*(nModes-1)/2;
  printmat(doublesVector2,numPsi,nPairs,numStates*numStates,1.0);
  printf("Triples\n");
  double *triplesVector2 = &integrals[1][0];
  int nTriples = (nModes)*(nModes-1)*(nModes-2)/6;
  printmat(triplesVector2,numPsi,nTriples,numStates*numStates*numStates,1.0);
  printf("Quadruples\n");
  double *quadruplesVector2 = &integrals[2][0];
  int nQuadruples = (nModes)*(nModes-1)*(nModes-2)*(nModes-3)/24;
  printmat(quadruplesVector2,numPsi,nQuadruples,numStates*numStates*numStates*numStates,1.0);

}

void Mp2Corr::fillCorrectionMatrices(int excitationLevel, int psiIndex, int potIndex, int tupleIndex, int startIndex, vector<int> diff, vector<double>& excitedEnergies) {
  if(pot[potIndex]->dim+1 >= excitationLevel) { //Note: this condition needs to change if we are using different degenerate subspaces from Gerber
    for(int i=startIndex ; i<nModes ; i++) {
      vector<int> diffCopy = diff; //copy diff
      diffCopy.push_back(i);

      for(int j=minState ; j<=maxState ; j++) { 
        dof[i]->setKet(j); 
        if(integralIsNonZero(diffCopy,tuples[potIndex][tupleIndex])) {
          int tupleStart = diffCopy.size()-excitationLevel;
          vector<int> indices(&diffCopy[tupleStart],&(*diffCopy.end())); 

          int index = getIndex(indices,psiIndex,excitationLevel);
          integrals[excitationLevel-2][index] += pot[potIndex]->integrateTuple(tupleIndex, false);

          //calculate denominator
          if(denominators[excitationLevel-2][index] == 0) { 
            double denominator = excitedEnergies[psiIndex+1]-excitedEnergies[0]; //En
            for(int a=0 ; a<indices.size() ; a++)
              denominator -= dof[indices[a]]->getEModal()-dof[indices[a]]->getEModal(0); //Em 
            denominators[excitationLevel-2][index] = denominator; 
          }
        }
        fillCorrectionMatrices(excitationLevel+1, psiIndex, potIndex, tupleIndex, i+1, diffCopy, excitedEnergies);
        dof[i]->setKet(0);
      }//state
    }//mode

  }//if 
}

void Mp2Corr::getSecondOrderCorr(vector<double>& corrections) {
  if(numPsi == 0) {
    std::cout << "No calculations have been done since the last clear.\n";
    return;
  } 
  double* coeff = new double[numPsi*numPsi];
  for(int i=0 ; i<numPsi ; i++) 
    for(int j=0 ; j<numPsi ; j++) 
      if(i==j)
        coeff[i*numPsi+j] = 1.0;
      else
        coeff[i*numPsi+j] = 0.0;
  getSecondOrderCorr(corrections,coeff); 
  delete[] coeff;
}

void Mp2Corr::getSecondOrderCorr(vector<double>& corrections, double* coeff) {
  //Check that calculateIntegrals was called before this and no clear
  if(numPsi == 0) {
    std::cout << "No calculations have been done since the last clear.\n";
    return;
  } 

  //VMP2 Corrections square the integrals and include denominator
  for(int i=0 ; i<numPsi ; i++) {
    int blockSize = singles.size()/numPsi;
    for(int j=0 ; j<blockSize ; j++) {
      double temp = 0.0;
      for(int k=0 ; k<numPsi ; k++) {
        temp += coeff[i*numPsi+k]*singles[k*blockSize+j];
      }
      if(temp != 0.0 && singlesDenom[i*blockSize+j] != 0.0)
        corrections[i] += temp*temp/singlesDenom[i*blockSize+j];
    }
    for(int j=0 ; j<integrals.size() ; j++) {
      blockSize = integrals[j].size()/numPsi;
      for(int k=0 ; k<blockSize ; k++) {
        double temp = 0.0;
        for(int l=0 ; l<numPsi ; l++) {
          temp += coeff[i*numPsi+l]*integrals[j][l*blockSize+k];
        }
        if(temp != 0.0 && denominators[j][i*blockSize+k] != 0.0)
          corrections[i] += temp*temp/denominators[j][i*blockSize+k];
      }
    } 
  }
}

void Mp2Corr::clear() {
  singles.clear();
  singlesDenom.clear();
  integrals.clear();
  denominators.clear();
  numPsi = 0;
}

//============Helper Methods==========
int Mp2Corr::getIndex(vector<int> indices, int psiIndex, int excitationLevel) {
  //excited bra block
  int index = psiIndex; 
  for(int i=0 ; i<excitationLevel ; i++) { 
    index *= nModes-i; 
    index /= i+1;
  }
  index *= pow(numStates,excitationLevel);

  //excited ket tuple
  index += tupleIndexDriver(indices,nModes)*pow(numStates,excitationLevel);

  //state combinations in tuple
  for(int i=0 ; i<indices.size() ; i++) {
    index += (dof[indices[i]]->getKet()-minState)*pow(numStates,excitationLevel-(i+1));
  }

  return index;
}

bool Mp2Corr::integralIsNonZero(vector<int> diff, vector<int> tuple) {
  for(int n=0 ; n<diff.size() ; n++) {
    if(dof[diff[n]]->getBra() != dof[diff[n]]->getKet()) {
      bool found = false;
      for(int n2=0 ; n2<tuple.size() ; n2++) {
        if(diff[n]==tuple[n2]) {
          found = true;
          break;
        }
      }
      if(!found) 
        return false;
    }          
  }
  return true;
}

//in cases of single excitations only, check if integral is 0. if a difference is only observed in a single mode, brillouin's theorem applies. (State is an eigenvector of the original state's Fock matrix)
//This is not related to Brillouin's: but if numDiff == 0, then the state is the same as the original state of interest, which should not be included. 
bool Mp2Corr::brillouin() {
  int numDiff = 0;
  for(int i=0 ; i<dof.size() ; i++) {
    if(dof[i]->getBra() != dof[i]->getKet())
      numDiff++; 
  } 
  return (numDiff == 1 || numDiff == 0);
}
