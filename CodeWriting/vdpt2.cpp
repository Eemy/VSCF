#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <stdbool.h>
#include "../UtilityFunc/aux.h"
#include "../Modules/Mode.h"
#include "../Modules/Mp2Corr.h"
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

void readin(std::vector<Mode*>& dof, std::vector<double>& freq, int N, int nPoints, int conv, int subspace);
bool checkConvergence(std::vector<Mode*> dof, double energy, int conv);
void print(FILE* script, std::string line);

double prevEnergy = 0.0;

int main(int argc, char* argv[]) {
  const int defaultLength = 10;
  int maxIter = 500;

  if(argc < defaultLength) { 
    printf("Error: <Nmodes> <Nquad> <1-Roothaan 2-Diis> <Diis_subspace> <EnergyFile> <DipoleXFile> <DipoleYFile> <DipoleZFile> <CouplingDegree> [<EnergyFile> <Dx> <Dy> <Dz> <dim> ...]\n");
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
  int diis_subspace = atoi(argv[4]);
  std::vector<std::string> potFileNames;
  potFileNames.push_back(argv[5]); //arg4
  std::vector<std::string> dipFileNames;
  dipFileNames.push_back(argv[6]); //arg5
  dipFileNames.push_back(argv[7]); //arg6
  dipFileNames.push_back(argv[8]); //arg7
  std::vector<int> potDims;
  potDims.push_back(atoi(argv[9])); //arg8

/////////////////////////Create Mode, EigSolver, Potential Objects////////////////////////
  std::vector<Mode*> dof;
  std::vector<double> freq;
  readin(dof,freq,nModes,nPoints,conv,diis_subspace); 
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
//    dof[i]->setHarmonic();
  }
//===================VMP2 Corrections to GS====================

  Mp2Corr correlationCalc(dof,pot,potIterators);
  std::vector<std::pair<std::vector<int>,std::vector<int>>> psi_m;

  std::vector<int> mode;
  std::vector<int> state(nModes,0);
  for(int i=0 ; i<nModes ; i++)
    mode.push_back(i);
  psi_m.push_back(std::make_pair(mode,state));//every mode at 0
  std::vector<double> energies;
  energies.push_back(excitedEnergies[0]); //GS energy
  energies.push_back(excitedEnergies[0]); //energy of state to be corrected (same state in this case)
  correlationCalc.calculateIntegrals(psi_m,energies);
  std::vector<double> mp2CorrGS(1);
  correlationCalc.getSecondOrderCorr(mp2CorrGS);
  
  printf("MP2 correction is (a.u) %15.12f\n",mp2CorrGS[0]);
  printf("MP2 correction is (a.u) %15.4f\n",mp2CorrGS[0]*au_to_wn);
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
  
  //Find single excitation energies (note: only works if maxQuanta=1 -> Gerber)
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
//need to save where the corresponding eigvec coefficients are for CI coeff and mp2 correction summing
//=========================End VCIS============================
/*  //read rmass for tests
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
*/
//======================VMP2 Corrections=======================
  correlationCalc.clear();
  psi_m.clear();
  for(int i=0 ; i<nModes ; i++) {
    std::vector<int> mode{i};
    std::vector<int> state{1};
    psi_m.push_back(std::make_pair(mode,state));
  }
  correlationCalc.calculateIntegrals(psi_m,excitedEnergies);
  std::vector<double> mp2Corr(nModes);
  correlationCalc.getSecondOrderCorr(mp2Corr,CI);

  printf("MP2 correction is (a.u) %15.12f\n",mp2Corr[0]);
  printf("MP2 correction is (a.u) %15.4f\n",mp2Corr[0]*au_to_wn);
  printf("MP2 correction is (a.u) %15.12f\n",mp2Corr[1]);
  printf("MP2 correction is (a.u) %15.4f\n",mp2Corr[1]*au_to_wn);
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
    fprintf(results," % -15.4f % -15.4f % \n", freq[i],(excitedEnergies[i+1]-excitedEnergies[0]+mp2Corr[i]+mp2CorrGS[0])*(au_to_wn));
//    fprintf(results," % -15.4f % -15.4f % \n", freq[i],(excitedEnergies[i+1]-excitedEnergies[0])*(au_to_wn));

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
void readin(std::vector<Mode*>& dof, std::vector<double>& freq, int N, int nPoints, int conv, int subspace) {
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
    dof.push_back(new Mode(freq[i],nPoints,conv,subspace));
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
