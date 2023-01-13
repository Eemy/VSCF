#include "Mode.h"
#include "Potential.h"

Mp2Corr::Mp2Corr(std::vector<Mode*> _dof, std::vector<Potential*> _pot, int maxDim) { 
  dof = _dof;
  pot = _pot;

  maxState = 4;
  minState = 1;
  numStates = maxState-minState+1;
  nModes = dof.size();

  singles.resize(nModes*nModes*numStates);
  singlesDenom.resize(nModes*nModes*numStates);

  for(int i=0 ; i<maxDim ; i++) {
    int excitationLevel = i+2;
    int size = nModes;
    for(int j=0 ; j<excitationLevel ; j++) { //computes unique tuples 
      size *= nModes-j; 
      size /= j+1;
    }
    size *= pow(numStates,excitationLevel);
    std::vector<double> integral(size);
    std::vector<double> denominator(size);
    integrals.push_back(integral);
    denominators.push_back(denominator);
  } 
}

Mp2Corr::calculateIntegrals() {
  
}

void fillCorrectionMatrices(Potential *pot, int minState, int maxState, std::vector<Mode*>& dof, int excitationLevel, std::vector<std::vector<double>>& integrals,  std::vector<int> tuple, int tupleIndex, int startIndex, std::vector<int> diff, std::vector<std::vector<double>>& denominators, std::vector<double>& excitedEnergies) {
  if(pot->dim+1 >= excitationLevel) {
    for(int i=startIndex ; i<dof.size() ; i++) {
      std::vector<int> diffCopy = diff; //copy diff
      diffCopy.push_back(i);
      for(int j=minState ; j<=maxState ; j++) { 
        dof[i]->setKet(j); 
        if(integralIsNonZero(diffCopy,tuple,dof)) {
          int index = getIndex(diffCopy,minState,maxState,excitationLevel,dof);
          integrals[excitationLevel-2][index] += pot->integrateTuple(tupleIndex, false);
          if(denominators[excitationLevel-2][index] == 0) { 
            std::vector<int> indices(&diffCopy[1],&(*diffCopy.end()));
            double denominator = excitedEnergies[diffCopy[0]+1]-excitedEnergies[0]; //En
            for(int a=0 ; a<indices.size() ; a++)
              denominator -= dof[indices[a]]->getEModal()-dof[indices[a]]->getEModal(0); //Em 
            denominators[excitationLevel-2][index] = denominator; 
          }
        }
        fillCorrectionMatrices(pot, minState, maxState, dof, excitationLevel+1, integrals, tuple, tupleIndex, i+1, diffCopy, denominators, excitedEnergies);
        dof[i]->setKet(0);
      }//state
    }//mode

  } 
}

int getIndex(std::vector<int> diff, int minState, int maxState, int excitationLevel, std::vector<Mode*> dof) {
  int numStates = maxState-minState+1;
  int nModes = dof.size();

  //excited bra block
  int index = diff[0]; //first element is the excited mode 
  for(int i=0 ; i<excitationLevel ; i++) { 
    index *= nModes-i; 
    index /= i+1;
  }
  index *= pow(numStates,excitationLevel);

  //excited ket tuple
  std::vector<int> indices(&diff[1],&(*diff.end()));
  index += tupleIndexDriver(indices,nModes)*pow(numStates,excitationLevel);

  //state combinations in tuple
  for(int i=0 ; i<indices.size() ; i++) {
    index += (dof[indices[i]]->getKet()-minState)*pow(numStates,excitationLevel-(i+1));
  }

  return index;
}

bool integralIsNonZero(std::vector<int> diff, std::vector<int> tuple, std::vector<Mode*>& dof) {
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

