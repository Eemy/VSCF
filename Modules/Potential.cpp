#include "Potential.h"
#include "Mode.h"
#include "cpot.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <unordered_map>
#include <string>
#include <vector>

std::string appendKey(int modeIndex) {
  if(modeIndex < 10)
    return "00"+std::to_string(modeIndex);
  if(modeIndex < 100)
    return "0"+std::to_string(modeIndex);
  return std::to_string(modeIndex);
}

//////////////////////////////////////////////START CLASS//////////////////////////////////////////////////
Potential::Potential(double* _pot, int _minDim, int _dimension, int _nModes, int _nPoints, int _totalLength, Mode** _modes) {
  dim = _dimension;
  if(_minDim == 1 && dim != 1) {
    minDim = 2;
  } else {
    minDim = _minDim;
  }
  fullPot = _pot;
  nModes = _nModes;
  nPoints = _nPoints;
  totalLength = _totalLength;
  potLength = (int) pow(_nPoints,_dimension);
  for(int i=0 ; i<dim ; i++) {
    subspaces.push_back((int) pow(nPoints,dim-i-1));
  }

//Set up Tuplets
  tuplets.reserve(totalLength/potLength);
  std::vector<int> iterations;
  for(int i=0 ; i<dim ; i++) {iterations.push_back(i);}
  int counter = 0;
  fillChecker(counter,0,iterations,_modes);
  for(int i = 0 ; i< totalLength/potLength ; i++) {
    tuplets[i]->setUpDriver(tupletChecker);
  }  
}

Potential::~Potential() {
  for(int i=0 ; i<totalLength/potLength ; i++) {
    delete tuplets[i];
  }
  delete[] fullPot;
}
//================================================================
//========================TUPLET NESTED CLASS=====================
//================================================================
Potential::Tuplet::Tuplet(int _minSubDim, int _dim, int _nPoints, std::vector<double> _pot, std::vector<Mode*> _modes, std::vector<int> _indices) {
  modeSubset = _modes;
  pot = _pot;
  indices = _indices;
  couplingRemover = new cpot(_nPoints);

  minSubDim = _minSubDim;
  dim = _dim;
  nPoints = _nPoints;
  potLength = (int) pow(nPoints,dim);
  for(int i=0 ; i<dim ; i++) {
    subspaces.push_back((int) pow(nPoints,dim-1-i));
  }
  if(dim > 2) {
   subTuplets.resize(dim-2);
  } 
}

Potential::Tuplet::~Tuplet() {
  for(int i=0 ; i<subTuplets.size() ; i++) {
    for(int j=0 ; j<subTuplets[i].size() ; j++) {
      delete subTuplets[i][j];
    }
  }
  delete couplingRemover;
}

//===================================================Nest within a nest=====================================================
Potential::Tuplet::SubTup::SubTup(std::vector<int> _subspaceSet, std::vector<Mode*> _modeSet, std::vector<double> _potSet) {
  subspaceSet = _subspaceSet;
  modeSet = _modeSet;
  potSet = _potSet; 
}
Potential::Tuplet::SubTup::~SubTup() {}
//=====================================================END================================================================

void Potential::Tuplet::setUpDriver(std::unordered_map<std::string,int>& checker) {
  setUp(dim,pot,checker,indices,subspaces,modeSubset);
  couplingRemover->get_coupling(pot, dim);
}

void Potential::Tuplet::setUp(int currDim, std::vector<double>& subPot, std::unordered_map<std::string,int>& checker, std::vector<int>& indicesInvolved, std::vector<int>& subspacesInvolved, std::vector<Mode*>& modes) {
  if(currDim > minSubDim) {
    for(int i=0 ; i<currDim; i++) {
      std::string key = "";
      std::vector<Mode*> modesInvolved;
      std::vector<int> newIndices;
      std::vector<int> newSubspaces;
      //create subtuplet identifiers
      for(int j=0 ; j<currDim; j++) {
        if(j>0)
          newSubspaces.push_back(subspacesInvolved[j]);
        if(j!=i) {
          key += appendKey(indicesInvolved[j]);
          modesInvolved.push_back(modes[j]);
          newIndices.push_back(indicesInvolved[j]);
        }
      }
      //check if the subtuplet has already been created under another tuplet 
      if(checker.count(key) == 0) {
        checker[key] = 1;
        std::vector<double> sliceNeeded(potLength/((int)pow(nPoints,dim-currDim+1)));
        getSlice(currDim,i,(nPoints-1)/2,subPot,sliceNeeded,subspacesInvolved);
        setUp(currDim-1,sliceNeeded,checker,newIndices,newSubspaces,modesInvolved);
        couplingRemover->get_coupling(sliceNeeded,currDim-1);
        subTuplets[currDim-3].push_back(new SubTup(newSubspaces,modesInvolved,sliceNeeded));
      }//if checker[key]
    }//end for i
  }
}

void Potential::Tuplet::getSlice(int dim, int modeSelection, int modeSlice, std::vector<double>& subPot, std::vector<double>& sliceNeeded, std::vector<int>& subspacesInvolved) {
  int halfPoint = (nPoints-1)/2;
  int center = ((int)pow(nPoints,dim)-1)/2;
  std::vector<int> otherSubspaces;

  //Find start of right slice
  int startIndex = center+(modeSlice-halfPoint)*subspacesInvolved[modeSelection];
  for(int i=0 ; i<dim ; i++) {
    if(i != modeSelection) {
      startIndex -= subspacesInvolved[i]*halfPoint;
      otherSubspaces.push_back(subspacesInvolved[i]);
    }
  }
  //Copy potential slice into a new array
  int iterations[20] = {0}; 
  fillPotential(dim-1, dim-1, sliceNeeded, subPot, startIndex, iterations, otherSubspaces);
}

void Potential::Tuplet::fillPotential(int loopsNeeded, int currDim, std::vector<double>& sliceNeeded, std::vector<double>& subPot, int startIndex, int* iters, std::vector<int>& otherSubSpace) {
  int recursionLevel = currDim-loopsNeeded;
  if(loopsNeeded == 1) {
    int indexStart = 0;
    int increment = 0;
    for(int i=0 ; i<currDim ; i++) {
      indexStart += (int)pow(nPoints,currDim-1-i)*iters[i];
      increment += otherSubSpace[i]*iters[i];
    }
    for(int i=0 ; i<nPoints; i++) {
      sliceNeeded[indexStart+i] = subPot[startIndex+increment+otherSubSpace[currDim-1]*i];
    }
  } else {
    for(int i=0 ; i<nPoints ; i++) {
      fillPotential(loopsNeeded-1, currDim, sliceNeeded, subPot, startIndex, iters, otherSubSpace);
      iters[recursionLevel]++;
      for(int j=recursionLevel+1 ; j<currDim ; j++)
        iters[j] = 0;
    }
  }
}

double Potential::Tuplet::getEffVIntegral(int modeSelection, int modeSlice) {
  double integralValue = 0.0;
  //Get the contribution from the Tuplet itself
  int weightMarkers[20] = {0};
  std::vector<double> tupletSlice(potLength/nPoints);
  std::vector<Mode*> modesNeeded;
  for(int i=0 ; i<modeSubset.size() ; i++) {
    if(i != modeSelection) 
      modesNeeded.push_back(modeSubset[i]);
  }
  getSlice(dim, modeSelection, modeSlice, pot, tupletSlice, subspaces);
  integralValue += getIntegral(0, tupletSlice.size(), dim-1, dim-1, modesNeeded, tupletSlice, weightMarkers);

  //include condition for 2mode VSCF
  if(dim > 2) {
    //loop through all subtuplets to check for selected mode
    Mode* modeSelected = modeSubset[modeSelection];
    for(int i=0 ; i<dim-2 ; i++) {
      for(int j=0 ; j<subTuplets[i].size() ; j++) {
        int indexOfMode = -1;
        for(int k=0 ; k<i+2 ; k++) {
          if(subTuplets[i][j]->modeSet[k] == modeSelected) {
            indexOfMode = k;
            break;
          }
        }
        //if tuplet contains the selected mode...integrate appropriate slice
        if(indexOfMode != -1) {
          int weightMarkers2[20] = {0};
          std::vector<double> sliceNeeded((int)pow(nPoints,i+1)); 
          std::vector<Mode*> modesInvolved;
          for(int k=0 ; k<i+2 ; k++) {
            if(k != indexOfMode)  
              modesInvolved.push_back(subTuplets[i][j]->modeSet[k]);
          }
          getSlice(i+2, indexOfMode, modeSlice, subTuplets[i][j]->potSet, sliceNeeded, subTuplets[i][j]->subspaceSet);
          integralValue += getIntegral(0,(int)pow(nPoints,i+1),i+1,i+1, modesInvolved, sliceNeeded, weightMarkers2);  
        }//if indexOfMode
      }//for j
    }//for i
  }//if dim>2
  return integralValue;
}

double Potential::Tuplet::getTotalIntegral(int excitedIndex) {
  if(excitedIndex != -1)
    modeSubset[excitedIndex]->setExcited(true);  

  double integralValue = 0.0;
  //obtain Tuplet's contribution
  int weightMarkers[20] = {0};
  integralValue += (dim-1)*getIntegral(0, potLength, dim, dim, modeSubset, pot, weightMarkers);
  
  //include condition for 2mode VSCF
  if(dim > 2) {
    //loop through all subtuplets
    for(int i=0 ; i<dim-2 ; i++) {
      for(int j=0 ; j<subTuplets[i].size() ; j++) {
        int weightMarkers2[20] = {0};
        integralValue += (i+1)*getIntegral(0,(int)pow(nPoints,i+2),i+2,i+2,subTuplets[i][j]->modeSet, subTuplets[i][j]->potSet,weightMarkers2); 
      }
    }
  }
  
  if(excitedIndex != -1)
    modeSubset[excitedIndex]->setExcited(false);
  return integralValue;
}

double Potential::Tuplet::getIntegral(int startIndex, int endIndex, int couplingDegree, int currentDim, std::vector<Mode*>& psi, std::vector<double>& potential, int* weightMarkers) {
  int potSize = endIndex - startIndex;
  double integralValue = 0.0;
  //Base Case: 1D
  if(potSize == nPoints) {
    for(int i=0 ; i<nPoints ; i++) {
      integralValue += psi[couplingDegree-1]->getIntegralComponent(i)*potential[startIndex+i]*psi[couplingDegree-1]->getWeight(i);
    }
    for(int i=0 ; i<couplingDegree-1 ; i++) {
      integralValue *= psi[i]->getIntegralComponent(weightMarkers[i])*psi[i]->getWeight(weightMarkers[i]);
    }
  } else {
    int newPotSize = potSize/nPoints;
    for(int i=0 ; i<nPoints ; i++) {
      int weightMarkerCopy[20] = {0};
      for(int j=0 ; j<couplingDegree ; j++) {
        weightMarkerCopy[j] = weightMarkers[j];
      }
      weightMarkerCopy[couplingDegree-currentDim] = i;
      integralValue += getIntegral(startIndex+i*newPotSize, startIndex+(i+1)*newPotSize, couplingDegree,
                      currentDim-1, psi, potential, weightMarkerCopy);
    }
  }
  return integralValue;
}
//================================================================================================
//==================================END TUPLET NESTED CLASS=======================================
//================================================================================================
double Potential::integralDriver(int index, int modeSelection, int modeSlice) {
  return tuplets[index]->getEffVIntegral(modeSelection,modeSlice);
}

double Potential::getVMinus() {
  double Vminus = 0.0;
  for(int i=0 ; i<totalLength/potLength ; i++) {
    Vminus += tuplets[i]->getTotalIntegral(-1);
  } 
  return Vminus;
}

double Potential::getDipole(int index, int excitedModeIndex) {
  return tuplets[index]->getTotalIntegral(excitedModeIndex);
}

double Potential::integrateSlice(Mode* mode, double* slice, bool excite) {
  if(excite)
    mode->setExcited(true);   

  int weightMarkers[20] = {0};
  std::vector<Mode*> psi;
  psi.push_back(mode);
  std::vector<double> potential;
  for(int i=0 ; i<nPoints ; i++)
    potential.push_back(slice[i]);    

  //It doesn't matter which tuple is used, just need integrating function
  double integralVal = tuplets[0]->getIntegral(0,nPoints,1,1,psi,potential,weightMarkers);  
  if(excite)
    mode->setExcited(false);
  return integralVal;
}
 
double** Potential::get1DSlices() {
  //this entire 2d array needs to be freed elsewhere. DON'T LET THE MEMORY LEAK HAPPEN
  double** slices = new double*[nModes];
  for(int i=0 ; i<nModes ; i++)
    slices[i] = new double[nPoints];

  int center = (potLength-1)/2;
  int halfPoint = (nPoints-1)/2;
  for(int i=0 ; i<nModes ; i++) {
    for(int j=0 ; j<nPoints ; j++) {
      int index = j-halfPoint;
      if(i<dim) {
        slices[i][j] = fullPot[center+index*subspaces[i]];
      } else {
        slices[i][j] = fullPot[(i-dim+1)*potLength+center+index*subspaces[dim-1]];
      }
    }
  }   
  return slices;
}

void Potential::fillChecker(int& counter, int recursionLevel, std::vector<int>& iter, Mode **_modes) {
  if(recursionLevel == dim-1) {
    while(iter[recursionLevel] < nModes) {
      std::vector<Mode*> modeSubset(dim);
      std::vector<int> indices(dim);
      for(int i=0 ; i<dim ; i++) {
        modeSubset[i] = _modes[iter[i]];
        indices[i] = iter[i];
      }
      std::vector<double> tupletPot(potLength);
      for(int i=0 ; i<potLength ; i++) {
        tupletPot[i] = fullPot[counter*potLength+i];     
      }
      tuplets.push_back(new Tuplet(minDim,dim,nPoints,tupletPot,modeSubset,indices));    
      iter[recursionLevel]++;
      counter++;
    }
//Recursive Cases
  } else { 
    while(iter[recursionLevel] < nModes-iter.size()+recursionLevel+1) {
      fillChecker(counter,recursionLevel+1,iter,_modes);
      iter[recursionLevel]++;
      for(int i=recursionLevel+1 ; i<dim ; i++) {
        iter[i] = iter[i-1]+1;     
      }
    }   
  }
}
