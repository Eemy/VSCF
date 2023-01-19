#include "Mode.h"
#include "Potential.h"
#include <vector>
#include <pair>
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
  //maxDim;
}

Mp2Corr::~Mp2Corr() {
}

Mp2Corr::calculateIntegrals(vector<pair<vector<int>,vector<int>>>psi_m, vector<double> excitedEnergies) {
  int numPsi = psi_m.size();  

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
  
  for(int i=0 ; i<numPsi ; i++) {
    for(int i2=0 ; i2<psi_m[i].first.size() ; i2++) 
      dof[psi_m[i].first[i2]]->setBra(psi_m[i].second[i2]);
    
    for(int j=0 ; j<tuples.size() ; j++) {
    for(int j2=0 ; j2<tuples[j].size() ; j2++) {

      //Single Excitations 
      for(int k=0 ; k<nModes ; k++) {
        std::vector<int> diff1;
        diff1.push_back(i); //change for psi_m
        diff1.push_back(k);
        for(int l=minState ; l<=maxState ; l++) {
          dof[k]->setKet(l);
          if(l>1 && integralIsNonZero(diff1,tuples[j][j2],dof)) { 
            double integralVal = pot[j]->integrateTuple(j2,false);
            int index = i*nModes*numStates+k*numStates+(l-minState);//change for psi_m
            //printf("Integral: %i %i State: %i, Tuple Num: %i Val: %.12f\n",i,k,l,j2,integralVal);   
            singles[index] += integralVal;
            if(singlesDenom[index] == 0.0)
              singlesDenom[index] = excitedEnergies[i+1]-(excitedEnergies[0]+(dof[k]->getEModal()-dof[k]->getEModal(0)));
          }

          //Double+ Excitations
          fillCorrectionMatrices(2, j, j2, k+1, diff1, excitedEnergies);
          dof[k]->setKet(0);
        }//l
      }//k
    }//j2
    }//j
    dof[i]->setBra(0); //change for psi_m
  }//i

  }    
}

Mp2Corr::clear() {
  singles.clear();
  singlesDenom.clear();
  integrals.clear();
  denominators.clear();
}

void Mp2Corr::fillCorrectionMatrices(int excitationLevel, int potIndex, int tupleIndex, int startIndex, std::vector<int> diff, std::vector<double>& excitedEnergies) {
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
        fillCorrectionMatrices(excitationLevel+1, potIndex, tupleIndex, i+1, diffCopy, excitedEnergies);
        dof[i]->setKet(0);
      }//state
    }//mode

  } 
}

int Mp2Corr::getIndex(std::vector<int> diff, int minState, int maxState, int excitationLevel, std::vector<Mode*> dof) {
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

bool Mp2Corr::integralIsNonZero(std::vector<int> diff, std::vector<int> tuple, std::vector<Mode*>& dof) {
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

