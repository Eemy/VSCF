#include "Potential.h"
#include "Mode.h"
#include "cpot.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <unordered_map>
#include <string>
#include <vector>
#include <fstream>
#include <stdbool.h>

std::string appendKey(int modeIndex) {
  if(modeIndex < 10)
    return "00"+std::to_string(modeIndex);
  if(modeIndex < 100)
    return "0"+std::to_string(modeIndex);
  return std::to_string(modeIndex);
}

//============================START CLASS============================
Potential::Potential(std::string fileName, int _minDim, int _dimension, int _nModes, int _nPoints) {
 dim = _dimension;
 if(_minDim == 1 && dim != 1) {
    minDim = 2;
  } else {
    minDim = _minDim;
  }
  nModes = _nModes;
  nPoints = _nPoints;
  potLength = (int) pow(_nPoints,_dimension);
  for(int i=0 ; i<dim ; i++)
    subspaces.push_back((int) pow(nPoints,dim-i-1));
  file = fileName;
//  readPot(_modes);
}

Potential::~Potential() {
  for(int i=0 ; i<tuplets.size() ; i++)
    delete tuplets[i];
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

//============================Nest within a nest============================
Potential::Tuplet::SubTup::SubTup(std::vector<int> _subspaceSet, std::vector<Mode*> _modeSet, std::vector<double> _potSet) {
  subspaceSet = _subspaceSet;
  modeSet = _modeSet;
  potSet = _potSet; 
}
Potential::Tuplet::SubTup::~SubTup() {}
//===================================END====================================

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

double Potential::Tuplet::getTotalIntegral(bool correction) {
  double integralValue = 0.0;
  //obtain Tuplet's contribution
  int weightMarkers[20] = {0};
  double temp = getIntegral(0, potLength, dim, dim, modeSubset, pot, weightMarkers);
  if(correction)
    temp *= (dim-1);
  integralValue += temp;
  
  //For 3D+ potentials: loop through all subtuplets
  for(int i=0 ; i<dim-2 ; i++) {
    for(int j=0 ; j<subTuplets[i].size() ; j++) {
      int weightMarkers2[20] = {0};
      double temp2 = getIntegral(0,(int)pow(nPoints,i+2),i+2,i+2,subTuplets[i][j]->modeSet, subTuplets[i][j]->potSet,weightMarkers2); 

      //Apply overlap integrals for modes not included 
      for(int k=0 ; k<dim ; k++) {
        bool found = false;
        for(int l=0 ; l<subTuplets[i][j]->modeSet.size() ; l++) {
          if(modeSubset[k] == subTuplets[i][j]->modeSet[l])
            found = true; 
        } 
        if(!found)
          temp2 *= modeSubset[k]->getOverlapEG();
      }

      if(correction)
        temp2 *= (i+1);
      integralValue += temp2; 
    }
  }
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
//============================================================================
//========================END TUPLET NESTED CLASS=============================
//============================================================================
double Potential::integralDriver(int index, int modeSelection, int modeSlice) {
  return tuplets[index]->getEffVIntegral(modeSelection,modeSlice);
}

double Potential::integrateTuple(int index, bool correction) {
  double integral = tuplets[index]->getTotalIntegral(correction);
  return integral;
}

double Potential::getDipole(int index, int excitedModeIndex) {
  return tuplets[index]->getTotalIntegral(excitedModeIndex);
}

double Potential::integrateSlice(Mode* mode, int modeIndex) {
  int weightMarkers[20] = {0};
  std::vector<Mode*> psi;
  psi.push_back(mode);

  //It doesn't matter which tuple is used, just need integrating function
  double integralVal = tuplets[0]->getIntegral(0,nPoints,1,1,psi,slices[modeIndex],weightMarkers);  
  return integralVal;
}
 
std::vector<std::vector<double>> Potential::get1DSlices() { return slices; }

//============================Big set-up function==========================
//Making the assumption that the lowest degree Potential is a complete set
std::vector<std::vector<int>> Potential::readPot(std::vector<Mode*>& dof, bool getSlices) {
  //Find out how many tuples are in the file
  //char tupleFile[20] = {0};
  char *tupleFile = new char[20];
  snprintf(tupleFile,19,"%i.dat",dim);
  //sprintf(tupleFile,"%i.dat",dim);
  //std::string tupleFile = std::to_string(dim) + ".dat";
  std::ifstream in(tupleFile,std::ios::in);
  if(!in) {
    printf("Error: %s could not be opened\n",tupleFile);
    exit(0);
  }
  delete[] tupleFile;
  std::string line;
  std::vector<std::vector<int>> tupleIndices;
  while(std::getline(in,line)) {
    //Parse tuple name to get modes and indices 
    std::vector<int> indices;
    std::string builder = "";
    for(auto &ch : line) {
      if(ch == '_') {
        indices.push_back(std::stoi(builder));
        builder = "";
      } else {
        builder += ch;
      } 
    }
    tupleIndices.push_back(indices);
  }
  in.close();
  
  //Read in actual potential, fill in tuple array
  in.open(file);
  double value = 0.0;
  int index = 0;
  int tupleIndex = 0;
  int potLength = pow(nPoints,dim); 
  std::vector<double> tuplePot(potLength);
  while(in >> value) {
    tuplePot[index%potLength] = value;
    index++;
    //Once potential is filled up...
    if(index%potLength == 0) {
      std::vector<Mode*> modes;
      for(int i=0 ; i<tupleIndices[tupleIndex].size() ; i++) 
        modes.push_back(dof[tupleIndices[tupleIndex][i]]);
      //Remove equilibrium energy
      double eqEnergy = tuplePot[(potLength-1)/2];
      for(int i=0 ; i<potLength ; i++) { 
        if(file.find("D") != 0) 
          tuplePot[i] -= eqEnergy;
        if(file.find("D") == 0)
          tuplePot[i] *= 0.393430307; //convert from debye_to_ea0 
      }

      tuplets.push_back(new Tuplet(minDim,dim,nPoints,tuplePot,modes,tupleIndices[tupleIndex]));  
      tupleIndex++;
    }
  }
  in.close();

 //Check potential lengths
  if(tupleIndex > tupleIndices.size() || index > tupleIndices.size()*potLength) {
    printf("%s is too long.\n",file.c_str());
    exit(0);
  } 
  if(tupleIndex < tupleIndices.size() || index < tupleIndices.size()*potLength) {
    printf("%s is too short.\n",file.c_str());
    exit(0);
  }

  //Prepare to read 1D Slices
  if(getSlices) {
    slices.resize(nModes); 
    for(int i=0 ; i<nModes ; i++)
      slices[i].resize(nPoints);

    int center = (potLength-1)/2;
    int halfPoint = (nPoints-1)/2;
    for(int i=0 ; i<nModes ; i++) {
      for(int j=0 ; j<nPoints ; j++) {
        int index = j-halfPoint;
        if(i<dim) {
          slices[i][j] = tuplets[0]->pot[center+index*subspaces[i]];
        } else {
          slices[i][j] = tuplets[(i-dim+1)]->pot[center+index*subspaces[dim-1]];
        }
    } }
  }   

  //Set up the tuples
  for(int i=0 ; i<tuplets.size() ; i++) 
    tuplets[i]->setUpDriver(tupletChecker);
  return tupleIndices;
} 
