#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdbool.h>
#include "../UtilityFunc/aux.h"
#include "../Modules/Mode.h"
#include "../Modules/Potential.h"
#include "../Modules/EigSolver.h"
///////////////////////////////CONSTANTS//////////////////////////////
#define pi       3.1415926535897932384626433832795
#define c_in_au  137.035999  //(a0/tau)
#define a0_to_cm 5.291772108e-9 
#define debye_to_ea0 0.393430307
#define Na       6.02214179E23 
#define au_to_wn 219474.6313708
//////////////////////////////////////////////////////////////////////

void fillCorrectionMatrices(Potential *pot, int minState, int maxState, std::vector<Mode*>& dof, int excitationLevel, std::vector<std::vector<double>>& integrals,  std::vector<int> tuple, int tupleIndex, int startIndex, std::vector<int> diff);
int getIndex(std::vector<int> diff, int minState, int maxState, int excitationLevel, std::vector<Mode*> dof);
bool integralIsNonZero(std::vector<int> diff, std::vector<int> tuple, std::vector<Mode*>& dof); 
void readin(std::vector<Mode*>& dof, std::vector<double>& freq, int N, int nPoints, int conv);
bool checkConvergence(std::vector<Mode*> dof, double energy, int conv);
void print(FILE* script, std::string line);

double prevEnergy = 0.0;

int main(int argc, char* argv[]) {
  const int defaultLength = 9;
  int maxIter = 500;

  if(argc < defaultLength) { 
    printf("Error: <Nmodes> <Nquad> <1-Roothaan 2-Diis> <EnergyFile> <DipoleXFile> <DipoleYFile> <DipoleZFile> <CouplingDegree> [<EnergyFile> <Dx> <Dy> <Dz> <dim> ...]\n");
    exit(0);
  }
  if((argc-defaultLength)%5 != 0) {
    printf("The number of args is invalid. Check your input and try again.\n");
    exit(0);
  }

  //SET-UP AND READ IN ARGS
  int nModes = atoi(argv[1]); //arg1
  int nPoints = atoi(argv[2]); //arg2
  int conv = atoi(argv[3]);
  if(conv != 1 && conv != 2) //default is roothaan for bad input
    conv = 1;
  std::vector<std::string> potFileNames;
  potFileNames.push_back(argv[4]); //arg4
  std::vector<std::string> dipFileNames;
  dipFileNames.push_back(argv[5]); //arg5
  dipFileNames.push_back(argv[6]); //arg6
  dipFileNames.push_back(argv[7]); //arg7
  std::vector<int> potDims;
  potDims.push_back(atoi(argv[8])); //arg8

/////////////////////////Create Mode, EigSolver, Potential Objects////////////////////////
  std::vector<Mode*> dof;
  std::vector<double> freq;
  readin(dof,freq,nModes,nPoints, conv); 
  EigSolver solver(nPoints,conv);

  std::vector<Potential*> pot;
  std::vector<Potential*> dip;
  std::vector<std::vector<std::vector<int>>> potIterators;
  std::vector<std::vector<std::vector<int>>> dipIterators; 

  //For first set of V and D
  pot.push_back(new Potential(potFileNames[0],1,potDims[0],nModes,nPoints));
  potIterators.push_back(pot[0]->readPot(dof,true));
  
  for(int i=0 ; i<3 ; i++) { 
    dip.push_back(new Potential(dipFileNames[i],1,potDims[0],nModes,nPoints));
    dipIterators.push_back(dip[i]->readPot(dof,true)); 
  }

  //For subsequent sets of V and D (5 represents number of args one set occupies in argv)
  if ((argc-defaultLength)%5 == 0 && (argc-defaultLength) > 0) {
    int len = (argc-defaultLength)/5;
    for(int i=1 ; i<=len ; i++) {
      potFileNames.push_back(argv[defaultLength+5*(i-1)]);
      potDims.push_back(atoi(argv[defaultLength+5*(i-1)+4]));
      pot.push_back(new Potential(potFileNames[i],potDims[i-1]+1,potDims[i],nModes,nPoints)); 
      potIterators.push_back(pot[i]->readPot(dof,false));
      for(int j=1 ; j<=3 ; j++) {
        dipFileNames.push_back(argv[defaultLength+5*(i-1)+j]);
        dip.push_back(new Potential(dipFileNames[i*3+(j-1)],potDims[i-1]+1,potDims[i],nModes,nPoints));
        dipIterators.push_back(dip[i*3+j-1]->readPot(dof,false));
      }
    } 
  } 
  //Get 1D slices
  std::vector<std::vector<double>> slices = pot[0]->get1DSlices();
  std::vector<double> excitedEnergies(nModes+1);
  std::vector<double> intensityComponents(3*nModes);
  std::vector<double> intensities(nModes);
  std::vector<double> overlaps(nModes);

  //Open results file once all set-up is completed
  FILE *results = fopen("eemVDPT2.dat","w");
//=========================Begin VSCF==============================
  //Prepare: eigensolver on pure 1D slices for each mode
  for(int i = 0 ; i< nModes ; i++) {
    prevEnergy += solver.solveMode(dof[i],slices[i],0,-1);//-1 prevents DIIS from occurring
  }
  //Compute effective potential integrals
  for(int iter = 0 ; iter< maxIter ; iter++) {
    std::vector<std::vector<double>> effV = pot[0]->get1DSlices();
    for(int i=0 ; i<potIterators.size() ; i++) {
      for(int j=0 ; j<potIterators[i].size() ; j++) {
        for(int k=0 ; k<potIterators[i][j].size() ; k++) {
          for(int l=0 ; l<nPoints ; l++) {
            int modeIndex = potIterators[i][j][k];
            effV[modeIndex][l] += pot[i]->integralDriver(j,k,l);
          }
        }//k loop: mode indices for tuple
      }//j loop: tuples 
    }//i loop: potentials

    //Compute VSCF Energy
    double energy = 0.0;
    for(int i = 0 ; i< nModes ; i++) {
      energy += solver.solveMode(dof[i],effV[i],0,iter);//iter instead of -1 allow DIIS to occur
    } 

    //Apply VSCF Energy Correction
    for(int i=0 ; i<pot.size() ; i++) {
      for(int j=0 ; j<potIterators[i].size() ; j++) {
        energy -= pot[i]->integrateTuple(j,true);
      }
    }

    //Check for Convergence
    if(checkConvergence(dof,energy,conv)) {
      fprintf(results,"Converged at iteration %d\n",iter+1);
      fprintf(results,"Ground-State VSCF Energy is: %.8f\n", energy*au_to_wn);
      excitedEnergies[0] = energy;
      break;
    } else {
      prevEnergy = energy;
    }
    if(iter == maxIter-1) {
      print(results,"Ground-State VSCF failed to converge.\n");
      excitedEnergies[0] = energy;
    }
  } 
//====================End Ground-State VSCF====================
  for(int i = 0 ; i< nModes ; i++) {
    dof[i]->setGroundState();
    dof[i]->setHarmonic();
  }
//===================VMP2 Corrections to GS====================
//====================End VMP2 Corrections=====================

//==================VCIS for all excited states================
  int maxQuanta = 1;
  double* CI = new double[nModes*nModes*maxQuanta*maxQuanta];
  for(int i=0 ; i<nModes*nModes*maxQuanta*maxQuanta ; i++) 
    CI[i] = 0.0;
  //Diagonal elements: VSCF + (En-E0) for the excited mode
  for(int i=0 ; i<maxQuanta ; i++) {
    for(int j=0 ; j<nModes ; j++) {
      CI[i*nModes*nModes*maxQuanta+j*nModes*maxQuanta+i*nModes+j] = excitedEnergies[0]+dof[j]->getEModal(i+1)-dof[j]->getEModal(0);
    }
  }

  //Off-diagonal elements
  for(int i=0 ; i<potIterators.size() ; i++) {
    for(int j=0 ; j<potIterators[i].size() ; j++) {
    //Iterate over each unique pair in tuple
      for(int k=0 ; k<potIterators[i][j].size() ; k++) {
        int firstMode = potIterators[i][j][k];
        for(int l=k+1 ; l<potIterators[i][j].size() ; l++) {
          int secondMode = potIterators[i][j][l];
          //Go through each single excitation block <100, <200, <300...
          for(int m=0 ; m<maxQuanta ; m++) {
            dof[firstMode]->setBra(m+1); 
            for(int n=0 ; n<maxQuanta ; n++) {
              dof[secondMode]->setKet(n+1);

              double integralVal = pot[i]->integrateTuple(j,false);
              CI[m*nModes*nModes*maxQuanta+firstMode*nModes*maxQuanta+n*nModes+secondMode] += integralVal; 
              CI[n*nModes*nModes*maxQuanta+secondMode*nModes*maxQuanta+m*nModes+firstMode] += integralVal; 

              dof[secondMode]->setKet(0); 
            }
            dof[firstMode]->setBra(0);
          }
        }//l loop: 2nd mode index for tuple
      }//k loop: mode indices for tuple
    }//j loop: tuples 
  }//i loop: potentials
 
  double* evals = new double[nModes*maxQuanta]; 
  diagonalize(CI,evals,nModes*maxQuanta);
  
  //Find single excitation energies
  int counter = 0;
  for(int i=0 ; i<nModes*maxQuanta ; i++) { 
    int maxIndex = 0;
    for(int j=0 ; j<nModes*maxQuanta ; j++) {
      if(fabs(CI[i*nModes*maxQuanta+maxIndex]) < fabs(CI[i*nModes*maxQuanta+j]))
        maxIndex = j; 
    }
    //The eigvec represents a single excitation
    if(maxIndex<nModes) {
      excitedEnergies[maxIndex+1] = evals[i];
      counter++;
    }
    if(counter==nModes)
      break;
  }
//=========================End VCIS============================
  //read rmass for tests
  std::ifstream in("rmass.dat",std::ios::in);
  std::vector<double> mass; 
 if(!in) {
    printf("Error: rmass.dat could not be opened\n");
    exit(0);
  }
  double val = 0.0;
  while(in >> val) {
    mass.push_back(val*1822.8884848961380701);
  }
//======================VMP2 Corrections=======================
  std::vector<double> mp2Corr(nModes);

  //perturbation is total - effV
  int maxState = 3;
  int minState = 1;
  int numStates = maxState-minState+1;

  //Perturbation integrals 
  std::vector<double> singles(nModes*nModes*numStates);
  double nPairs = nModes*(nModes-1)/2;
  std::vector<double> doubles(nModes*nPairs*numStates*numStates);
  double nTriples = nModes*(nModes-1)*(nModes-2)/6;
  std::vector<double> triples(nModes*nTriples*numStates*numStates*numStates); 
  for(int i=0 ; i<nModes ; i++) {
    dof[i]->setBra(1);
    for(int j=0 ; j<potIterators.size() ; j++) {
    for(int j2=0 ; j2<potIterators[j].size() ; j2++) {
      //Single Excitations 
      for(int k=0 ; k<nModes ; k++) {
        std::vector<int> diff1;
        diff1.push_back(i);
        diff1.push_back(k);
        //int singlesMin = minState>1 ? minState:2;//remove degenerate subspace
        for(int l=minState ; l<=maxState ; l++) {
          dof[k]->setKet(l);
          if(integralIsNonZero(diff1,potIterators[j][j2],dof)) { 
            double integralVal = pot[j]->integrateTuple(j2,false);
//            printf("Integral: %i %i State: %i, Tuple Num: %i Val: %.12f\n",i,k,l,j2,integralVal);   
            singles[i*nModes*numStates+k*numStates+(l-minState)] += integralVal;
          }

      //Double Excitations
      for(int k2=k+1 ; k2<nModes ; k2++) {
        std::vector<int> diff2;
        diff2.push_back(i);
        diff2.push_back(k);
        diff2.push_back(k2);
        for(int l2=minState ; l2<=maxState ; l2++) {
          dof[k2]->setKet(l2);
          if(integralIsNonZero(diff2,potIterators[j][j2],dof)) { 
//            printf("Integral: %i %i %i, Tuple Num: %i\n",i,k,l,j2);        
            double integralVal = pot[j]->integrateTuple(j2,false);
            auto start = diff2.begin()+1;
            auto end = diff2.end();
            std::vector<int> indices(diff2.size()-1);
            std::copy(start,end,indices.begin());
            doubles[i*nPairs*numStates*numStates+
                    tupleIndexDriver(indices,nModes)*numStates*numStates+
                    (l-minState)*numStates+
                    (l2-minState)] += integralVal;
          }


      //Triple Excitations
      for(int k3=k2+1 ; k3<nModes ; k3++) {
        std::vector<int> diff3;
        diff3.push_back(i);
        diff3.push_back(k);
        diff3.push_back(k2);
        diff3.push_back(k3);
        for(int l3=minState ; l3<=maxState ; l3++) {
          dof[k3]->setKet(l3);
          if(integralIsNonZero(diff3,potIterators[j][j2],dof)) { 
//            printf("Integral: %i %i %i, Tuple Num: %i\n",i,k,l,j2);        
            double integralVal = pot[j]->integrateTuple(j2,false);
            auto start = diff3.begin()+1;
            auto end = diff3.end();
            std::vector<int> indices(diff3.size()-1);
            std::copy(start,end,indices.begin());
            triples[i*nTriples*numStates*numStates*numStates+
                    tupleIndexDriver(indices,nModes)*numStates*numStates*numStates+
                    (l-minState)*numStates*numStates+
                    (l2-minState)*numStates+
                    (l3-minState)] += integralVal;
          }
          dof[k3]->setKet(0);
        }//l3
      }//k3

          dof[k2]->setKet(0);
        }//l2
      }//k2

          dof[k]->setKet(0);
        }//l
      }//k

    }//j2
    }//j
    dof[i]->setBra(0);
  }//i

  printf("\n");

  //Perturbation integrals 
  int maxDim = 0;
  for(int i=0 ; i<potDims.size() ; i++) {
    if(potDims[i] > maxDim)
      maxDim = potDims[i];
  }  

  std::vector<double> singles2(nModes*nModes*numStates);
  std::vector<std::vector<double>> integrals;
  for(int i=0 ; i<maxDim ; i++) {
    int excitationLevel = i+2;
    int size = nModes;
    for(int j=0 ; j<excitationLevel ; j++) { //computes unique tuples 
      size *= nModes-j; 
      size /= j+1;
    }
    size *= pow(numStates,excitationLevel);
    std::vector<double> integral(size);
    integrals.push_back(integral);
  } 
  for(int i=0 ; i<nModes ; i++) {
    dof[i]->setBra(1);
    for(int j=0 ; j<potIterators.size() ; j++) {
    for(int j2=0 ; j2<potIterators[j].size() ; j2++) {

      //Single Excitations 
      for(int k=0 ; k<nModes ; k++) {
        std::vector<int> diff1;
        diff1.push_back(i);
        diff1.push_back(k);
//        int singlesMin = minState>1 ? minState:2;//remove degenerate subspace
        for(int l=minState ; l<=maxState ; l++) {
          dof[k]->setKet(l);
          if(integralIsNonZero(diff1,potIterators[j][j2],dof)) { 
            double integralVal = pot[j]->integrateTuple(j2,false);
//            printf("Integral: %i %i State: %i, Tuple Num: %i Val: %.12f\n",i,k,l,j2,integralVal);   
            singles2[i*nModes*numStates+k*numStates+(l-minState)] += integralVal;
          }

          //Double+ Excitations
          fillCorrectionMatrices(pot[j], minState, maxState, dof, 2, integrals, potIterators[j][j2], j2, k+1, diff1);
          dof[k]->setKet(0);
        }//l
      }//k
    }//j2
    }//j
    dof[i]->setBra(0);
  }//i
  double mode1 = 1/sqrt(2*mass[0]*dof[0]->getOmega());
  double mode2 = 1/sqrt(2*mass[1]*dof[1]->getOmega());
  double mode3 = 1/sqrt(2*mass[2]*dof[2]->getOmega());
  double mode4 = 1/sqrt(2*mass[3]*dof[3]->getOmega());

  printf("Predicted: 1000|1122 %.12f\n",mode2*mode3*mode3*mode4*mode4*2);
  printf("Predicted: 0100|1122 %.12f\n",mode1*mode3*mode3*mode4*mode4*2);
  printf("Predicted: 0010|1212 %.12f\n",mode1*mode2*mode2*mode4*mode4*2);
  printf("Predicted: 0004|1221 %.12f\n",mode1*mode2*mode2*mode3*mode3*2);
/*
  printf("Singles\n");
  double *singlesVector = &singles[0];
  printmat(singlesVector,nModes,nModes,numStates,1.0);
  double *singlesVector2 = &singles2[0];
  printmat(singlesVector2,nModes,nModes,numStates,1.0);

  printf("Doubles\n");
  double *doublesVector = &doubles[0];
  printmat(doublesVector,nModes,nPairs,numStates*numStates,1.0);
  double *doublesVector2 = &integrals[0][0];
  printmat(doublesVector2,nModes,nPairs,numStates*numStates,1.0);

  printf("Triples\n");
  double *triplesVector = &triples[0];
  printmat(triplesVector,nModes,nTriples,numStates*numStates*numStates,1.0);
  double *triplesVector2 = &integrals[1][0];
  printmat(triplesVector2,nModes,nTriples,numStates*numStates*numStates,1.0);
*/
  
  
  printf("Quadruples\n");
  double *quadruplesVector2 = &integrals[2][0];
  int nQuadruples = (nModes)*(nModes-1)*(nModes-2)*(nModes-3)/24;
  printmat(quadruplesVector2,nModes,nQuadruples,numStates*numStates*numStates*numStates,1.0);
//===================End VMP2 Corrections======================

  //Print out all the transition frequencies
  print(results,"************************************************************************\n");
  print(results," VSCF ground-state energy: (cm^-1) \n");
  fprintf(results," % -15.4f \n", excitedEnergies[0]*(au_to_wn));
  print(results," \n");
  print(results," Transitions: (cm^-1) \n");
  print(results," \n");
  print(results,"  Harmonic        CIS(1)            Intensity(km/mol)\n");
  for(int i=0; i<nModes ; i++) {
//  fprintf(results," % -15.4f % -15.4f % -15.4f \n", freq[i],(excitedEnergies[i+1]-excitedEnergies[0])*(au_to_wn),intensities[i]);
    fprintf(results," % -15.4f % -15.4f % \n", freq[i],(excitedEnergies[i+1]-excitedEnergies[0])*(au_to_wn));

  }

  //DEALLOCATE
  for(int i=0 ; i<nModes ; i++) delete dof[i];
  for(int i=0 ; i<dip.size() ; i++) delete dip[i];
  for(int i=0 ; i<pot.size() ; i++) delete pot[i];
  delete[] evals;
  delete[] CI;

  return 0;
}//end main

//==========================HELPER METHODS=============================
void fillCorrectionMatrices(Potential *pot, int minState, int maxState, std::vector<Mode*>& dof, int excitationLevel, std::vector<std::vector<double>>& integrals,  std::vector<int> tuple, int tupleIndex, int startIndex, std::vector<int> diff) {
  if(pot->dim+1 >= excitationLevel) {
    for(int i=startIndex ; i<dof.size() ; i++) {
      std::vector<int> diffCopy = diff; //copy diff
      diffCopy.push_back(i);
      for(int j=minState ; j<=maxState ; j++) { 
        dof[i]->setKet(j); 
        if(integralIsNonZero(diffCopy,tuple,dof)) {
          if(excitationLevel == 3) {
            printf("%i %i %i %i |%i %i %i %i, Index: %i, Tuple: %i\n",dof[0]->getBra(),dof[1]->getBra(),dof[2]->getBra(),dof[3]->getBra(),dof[0]->getKet(),dof[1]->getKet(),dof[2]->getKet(),dof[3]->getKet(),getIndex(diffCopy,minState,maxState,excitationLevel,dof),tupleIndex);
          }
          integrals[excitationLevel-2][getIndex(diffCopy,minState,maxState,excitationLevel,dof)] += pot->integrateTuple(tupleIndex, false);
        }
        fillCorrectionMatrices(pot, minState, maxState, dof, excitationLevel+1, integrals, tuple, tupleIndex, i+1, diffCopy);
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

void readin(std::vector<Mode*>& dof, std::vector<double>& freq, int N, int nPoints, int conv) {
  //read in frequencies 
  std::ifstream in("freq.dat",std::ios::in);
  if(!in) {
    printf("Error: freq.dat could not be opened\n");
    exit(0);
  }
  double val = 0.0;
  while(in >> val) {
    freq.push_back(val);
  }

  //Check size of freq file
  if(freq.size() != N) {
    printf("freq.dat is the wrong size.\n");
    exit(0);
  }

  //Create Mode objects and return
  for(int i=0 ; i<N ; i++) {
    dof.push_back(new Mode(freq[i],nPoints,conv));
  }    
} 

bool checkConvergence(std::vector<Mode*> dof, double energy, int conv) {
  //Roothaan
  if(conv==1) {
    double diff = 0.0;
    for(int i=0 ; i<dof.size() ; i++) {
      double temp = dof[i]->computeMaxDiff();
      if(temp > diff)
        diff = temp;
    }
    return (diff < 1.0E-5) && (fabs(energy-prevEnergy)*au_to_wn <0.5);
  }
  //DIIS
  if(conv==2) {
    double max = 0.0;
    for(int i=0 ; i<dof.size() ; i++) {
      if(dof[i]->getDIISError() > max)
        max = dof[i]->getDIISError();
    }
    printf("ConvCheck: %.12f\n",max);
    return (max < 1.0e-10);
  }
}
void print(FILE* script, std::string line) {
  fprintf(script,line.c_str());
}
