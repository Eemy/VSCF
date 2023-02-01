#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdbool.h>
#include "../Modules/Mode.h"
#include "../Modules/Potential.h"
#include "../Modules/EigSolver.h"
///////////////////////////////CONSTANTS//////////////////////////////
#define pi       3.1415926535897932384626433832795
#define c_in_au  137.035999  //(a0/tau)
#define a0_to_cm 5.291772108e-9 
#define debye_to_ea0 0.393430307
#define Na       6.02214179E23 
//////////////////////////////////////////////////////////////////////
void readin(std::vector<Mode*>& dof, std::vector<double>& freq, int N, int nPoints, int conv);
bool checkConvergence(std::vector<Mode*> dof, double energy, int conv);
void print(FILE* script, std::string line);

double prevEnergy = 0.0;

int main(int argc, char* argv[]) {
  const int defaultLength = 9;
  int maxIter = 500;

  if(argc < defaultLength) { 
    printf("Error: <Nmodes> <Nquad> <1-Roothaan 2-Diis 3-GradDescent> <EnergyFile> <DipoleXFile> <DipoleYFile> <DipoleZFile> <CouplingDegree> [<EnergyFile> <Dx> <Dy> <Dz> <dim> ...]\n");
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
  if(conv != 1 && conv != 2 && conv != 3) //default is roothaan for bad input
    conv = 1;
  std::vector<std::string> potFileNames;
  potFileNames.push_back(argv[4]); //arg4
  std::vector<std::string> dipFileNames;
  dipFileNames.push_back(argv[5]); //arg5
  dipFileNames.push_back(argv[6]); //arg6
  dipFileNames.push_back(argv[7]); //arg7
  std::vector<int> potDims;
  potDims.push_back(atoi(argv[8])); //arg8

/////////////////Create Mode, EigSolver, Potential Objects///////////////////
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
  std::string fileName = "";
  if(conv==1)
     fileName = "eemVSCF_roothaan.dat";
  else if(conv==2)
     fileName = "eemVSCF_diis.dat";
  else
    fileName = "eemVSCF_gradDescent.dat";
  FILE *results = fopen(fileName.c_str(),"w");
//===============================Begin VSCF==================================
  //Prepare: eigensolver on pure 1D slices for each mode
  for(int i = 0 ; i< nModes ; i++) {
    prevEnergy += solver.solveMode(dof[i],slices[i],0,-1);//-1 prevents DIIS
    if(conv==2)
      dof[i]->saveCurrentDensity(-1);
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
      printf("Mode %i\n",i);
      if(iter==3)
        printf("checkpoint\n");
      energy += solver.solveMode(dof[i],effV[i],0,iter);//iter instead of -1 allow DIIS to occur
    } 

    //Apply VSCF Energy Correction
    for(int i=0 ; i<pot.size() ; i++) {
      for(int j=0 ; j<potIterators[i].size() ; j++) {
        energy -= pot[i]->integrateTuple(j,true);
      }
    }

    //CALL DIIS: solve for coefficients and obtain extrapolated density matrix
    if(conv==2) {
      solver.diis(dof,iter);
      //for(int i=0 ; i<dof.size() ; i++) dof[i]->saveCurrentDensity(iter);
    }

    //Check for Convergence
    if(checkConvergence(dof,energy,conv)) {
      fprintf(results,"Converged at iteration %d\n",iter+1);
      fprintf(results,"Ground-State VSCF Energy is: %.8f\n", energy*219474.6313708);
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
///////End Ground-State VSCF/////////
  for(int i = 0 ; i< nModes ; i++) {
    dof[i]->setGroundState();
    if(conv==2)
      dof[i]->resetSubspace();
  }
/////////Excited-State VSCF//////////
for(int z = 0 ; z< nModes ; z++) {
  //Prepare 1D slices (effV guess)
  prevEnergy = 0.0;
  dof[z]->setBra(1); //excite mode
  dof[z]->setKet(1);
  for(int i = 0 ; i< nModes ; i++) {
    if(i==z) {
      prevEnergy += solver.solveMode(dof[i],slices[i],1,-1);
    } else {
      prevEnergy += solver.solveMode(dof[i],slices[i],0,-1);
    }
    if(conv==2)
      dof[i]->saveCurrentDensity(-1);
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
      if(i==z) {
        energy += solver.solveMode(dof[i],effV[i],1,iter);
      } else {
        energy += solver.solveMode(dof[i],effV[i],0,iter);
      }
    }

    //Apply VSCF Energy Correction
    for(int i=0 ; i<pot.size() ; i++) {
      for(int j=0 ; j<potIterators[i].size() ; j++) {
        energy -= pot[i]->integrateTuple(j,true);
      }
    }
    if(z==1)
      printf("checkpoint\n");
    //CALL DIIS: solve for coefficients and obtain extrapolated density matrix
    if(conv==2) {
      solver.diis(dof,iter);
      //for(int i=0 ; i<dof.size() ; i++) dof[i]->saveCurrentDensity(iter);
    }

    //Check for Convergence
    if(checkConvergence(dof,energy,conv)) {
      fprintf(results,"Converged at iteration %d\n",iter+1);
      fprintf(results,"Mode %i Excited-State VSCF Energy is: %.8f\n",z,energy*219474.6313708);
      excitedEnergies[z+1] = energy;
      break;
    } else {
      prevEnergy = energy;
    }
    if(iter == maxIter-1) {
      fprintf(results,"Mode %i VSCF failed to converge.\n",z);
      excitedEnergies[z+1] = energy;
    }
  }
  dof[z]->setExcitedState();
  dof[z]->setBra(0);
  dof[z]->setKet(0);
  if(conv==2) {
    for(int i=0 ; i<dof.size() ; i++)
      dof[i]->resetSubspace(); 
  }
}//z loop: excited mode
////////End Excited-State VSCF///////

for(int i = 0 ; i< nModes ; i++) {
  dof[i]->useVSCFStates(true);
  overlaps[i] = dof[i]->getOverlapEG();
}

//DIPOLE CALCULATIONS
for(int comp = 0 ; comp< 3 ; comp++) {
  //Integrate all 1D pieces
  for(int a = 0 ; a< nModes ; a++) {
    dof[a]->setBra(1);
    intensityComponents[3*a+comp] += dip[comp]->integrateSlice(dof[a],a);
    dof[a]->setBra(0);
    for(int b = a+1 ; b< nModes ; b++) {
      intensityComponents[3*a+comp] += dip[comp]->integrateSlice(dof[b],b)*overlaps[a];
      intensityComponents[3*b+comp] += dip[comp]->integrateSlice(dof[a],a)*overlaps[b];
    }
  }

  for(int i=0 ; i<dipIterators.size()/3 ; i++) {
    for(int j=0 ; j<dipIterators[3*i+comp].size() ; j++) {
      for(int k=0 ; k<dipIterators[3*i+comp][j].size() ; k++) {
        int modeIndex = dipIterators[3*i+comp][j][k];
        dof[modeIndex]->setBra(1);
        dof[modeIndex]->setKet(0);
        intensityComponents[3*modeIndex+comp] += dip[3*i+comp]->getDipole(j); 
        dof[modeIndex]->setBra(0);
        dof[modeIndex]->setKet(0);
      }//k loop: mode indices for tuple
    }//j loop: tuples 
  }//i loop: potentials
}//comp loop: 0=x, 1=y, 2=z

for(int i=0 ; i<nModes ; i++) {
  for(int j=0 ; j<3 ; j++) {
    intensities[i] += intensityComponents[3*i+j]*intensityComponents[3*i+j]; 
  }
  intensities[i] *= (excitedEnergies[i+1]-excitedEnergies[0])*2.0*pi*Na/(3.0*c_in_au*c_in_au)*a0_to_cm*1.0E-5; //to km/mol
}

  //Print out all the transition frequencies
  print(results,"************************************************************************\n");
  print(results," VSCF ground-state energy: (cm^-1) \n");
  fprintf(results," % -15.4f \n", excitedEnergies[0]*(219474.6313708));
  print(results," \n");
  print(results," Transitions: (cm^-1) \n");
  print(results," \n");
  print(results,"  Harmonic        VSCF            Intensity(km/mol)\n");
  for(int i=0; i<nModes ; i++) {
  fprintf(results," % -15.4f % -15.4f % -15.4f \n", freq[i],(excitedEnergies[i+1]-excitedEnergies[0])*(219474.6313708),intensities[i]);
  }

  //DEALLOCATE
  for(int i=0 ; i<nModes ; i++) delete dof[i];
  for(int i=0 ; i<dip.size() ; i++) delete dip[i];
  for(int i=0 ; i<pot.size() ; i++) delete pot[i];

  return 0;
}

//==========================HELPER METHODS=============================
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
  if(conv==1 || conv==3) {
    double diff = 0.0;
    for(int i=0 ; i<dof.size() ; i++) {
      double temp = dof[i]->computeMaxDiff();
      if(temp > diff)
        diff = temp;
    }
    printf("Energy: %.12f\n",energy);
    return (diff < 1.0E-5) && (fabs(energy-prevEnergy)*219474.6313708 <0.5);
  }
  //DIIS
  if(conv==2) {
    double max = 0.0;
    for(int i=0 ; i<dof.size() ; i++) {
      if(dof[i]->getDIISError() > max)
        max = dof[i]->getDIISError();
    }
    printf("ConvCheck: %.12f\n",max);
    printf("Energy: %.12f\n",energy);
    return (max < 1.0e-14);
  }
}
void print(FILE* script, std::string line) {
  fprintf(script,line.c_str());
}
