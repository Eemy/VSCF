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

void readin(double* freq, double* mass, int N);
void readPot(std::string, double* potential, int length, int potLength);
bool checkConvergence(Mode** dof, double energy, int nModes);
int fact(int n);
void print(FILE* script, std::string line); 
int lengthCheck(int nMode, int coupling);

double prevEnergy = 0.0;

int main(int argc, char* argv[]) {
  if(argc!=6) printf("Error: <Nmodes> <Nquad> <couplingDegree> <#ofMethodsUsed> <damping>\n");
  else {
  //SET-UP
    int nModes = atoi(argv[1]);
    int nPoints = atoi(argv[2]);
    int couplingDegree = atoi(argv[3]);
    int iter = 100;
    int numMethods = atoi(argv[4]);    
    double b = atof(argv[5]);

    std::vector<int> potDims;
    std::vector<int> expectedLengths;
    int expectedLength = 0;
    std::vector<std::string> potFileNames;
    std::vector<double*> potentials;
    double *V;
//================Get info from user if more than one method is used====================
    if(numMethods > 1) {
      potDims.resize(numMethods);
      potFileNames.resize(numMethods);
      for(int i=0 ; i<numMethods ; i++) {
        std::cout << "File name of potential:";
        std::cin >> potFileNames[i];
        std::cout << "Dimension of potential: (should be in ascending order)";
        std::cin >> potDims[i];
      }
      char finalCheck = 'n';
      std::cout << "Last Q: Is the first potential the one you want to take 1D slices from?(y/n)";
      std::cin >> finalCheck;
      if(finalCheck != 'y') { 
        std::cout << "Go back and change it then!\n";
        exit(0);       
      }

      //Read in the potentials
      expectedLengths.resize(numMethods);
      potentials.resize(numMethods);
      for(int i=0 ; i<numMethods ; i++) {
        expectedLengths[i] = lengthCheck(nModes,potDims[i])*(int) pow(nPoints,potDims[i]);
        potentials[i] = new double[expectedLengths[i]];
        readPot(potFileNames[i],potentials[i],(int) pow(nPoints,potDims[i]),expectedLengths[i]);
      }
    } else { //if just a single method is used
      expectedLength = lengthCheck(nModes,couplingDegree)*(int)pow(nPoints,couplingDegree);
      V = new double[expectedLength];
      readPot("V.dat",V,(int) pow(nPoints,couplingDegree),expectedLength);
    }
//=====================================================================================

    EigSolver solver(nPoints);
    double* freq = new double[nModes];
    double* mass = new double[nModes];
    readin(freq, mass, nModes); 
    Mode** dof = new Mode*[nModes];
    for(int i=0 ; i<nModes ; i++) {
      dof[i] = new Mode(freq[i],mass[i],nPoints);
    }    
    double** effV = new double*[nModes];
    for(int i=0 ; i<nModes ; i++) {
      effV[i] = new double[nPoints];
    }
    double* excitedEnergies = new double[nModes+1];
    FILE *results = fopen("eemVSCF.dat","w");
//=====================================================================================
Potential pot(V,1,couplingDegree,nModes,nPoints,expectedLength,dof);
double** slices = pot.get1DSlices();

double** oldEffV = new double*[nModes];
for(int i=0 ; i<nModes ; i++) {
  oldEffV[i] = new double[nPoints];
}

//Prepare: eigensolver on pure 1D slices for each mode
for(int i = 0 ; i< nModes ; i++) {
prevEnergy += solver.solveMode(dof[i],slices[i],0);
}
for(int iter = 1 ; iter< 100 ; iter++) {
int counter = 0;
for(int i = 0 ; i< nModes ; i++) {
for(int j = 0 ; j< nPoints ; j++) {
effV[i][j] = slices[i][j];
}
}
for(int a = 0 ; a< nModes ; a++) {
for(int b = a+1 ; b< nModes ; b++) {
for(int c = b+1 ; c< nModes ; c++) {
for(int d = 0 ; d< nPoints ; d++) {
effV[a][d] += pot.integralDriver(counter, 0, d);
effV[b][d] += pot.integralDriver(counter, 1, d);
effV[c][d] += pot.integralDriver(counter, 2, d);
}
counter++;
}
}
}

/////Damping convergence:////
if(iter > 1) {
  for(int i=0 ; i<nModes ; i++) {
    for(int j=0 ; j<nPoints ; j++) {
      double newVal = effV[i][j];
      effV[i][j] = b*oldEffV[i][j]+(1-b)*newVal; 
    }
  }
}

for(int i=0 ; i<nModes ; i++) {
  for(int j=0 ; j<nPoints ; j++) {
    oldEffV[i][j] = effV[i][j];
  }
}
/////////////////////////////

double energy = 0.0;
for(int i = 0 ; i< nModes ; i++) {
energy += solver.solveMode(dof[i],effV[i],0);
}
energy -= pot.getVMinus();
if(checkConvergence(dof,energy,nModes)) {
fprintf(results,"Converged at iteration %d\n",iter);
fprintf(results,"Ground-State VSCF Energy is: %.8f\n", energy*219474.6313708);
excitedEnergies[0] = energy;
break;
} else {
prevEnergy = energy;
}
if(iter == 99) {
print(results,"VSCF failed to converge.\n");
excitedEnergies[0] = energy;
}}
/////////////////End Ground-State VSCF/////////////////////

//Save ground state for dipole calculations
for(int i=0 ; i<nModes ; i++) {
 dof[i]->setGroundState(); 
}

////////////////////Excited-State VSCF/////////////////////
for(int z = 0 ; z< nModes ; z++) {
prevEnergy = 0.0;
for(int i = 0 ; i< nModes ; i++) {
if(i==z) {
prevEnergy += solver.solveMode(dof[i],slices[i],1);
} else {
 prevEnergy += solver.solveMode(dof[i],slices[i],0);
}
}
for(int iter = 1 ; iter< 100 ; iter++) {
int counter = 0;
for(int i = 0 ; i< nModes ; i++) {
for(int j = 0 ; j< nPoints ; j++) {
effV[i][j] = slices[i][j];
}
}
for(int a = 0 ; a< nModes ; a++) {
for(int b = a+1 ; b< nModes ; b++) {
for(int c = b+1 ; c< nModes ; c++) {
for(int d = 0 ; d< nPoints ; d++) {
effV[a][d] += pot.integralDriver(counter, 0, d);
effV[b][d] += pot.integralDriver(counter, 1, d);
effV[c][d] += pot.integralDriver(counter, 2, d);
}
counter++;
}
}
}

/////Damping convergence:////
if(iter > 1) {
  for(int i=0 ; i<nModes ; i++) {
    for(int j=0 ; j<nPoints ; j++) {
      double newVal = effV[i][j];
      effV[i][j] = b*oldEffV[i][j]+(1-b)*newVal; 
    }
  }
}

for(int i=0 ; i<nModes ; i++) {
  for(int j=0 ; j<nPoints ; j++) {
    oldEffV[i][j] = effV[i][j];
  }
}
/////////////////////////////

double energy = 0.0;
for(int i = 0 ; i< nModes ; i++) {
if(i==z) {
energy += solver.solveMode(dof[i],effV[i],1);
} else {
energy += solver.solveMode(dof[i],effV[i],0);
}
}
energy -= pot.getVMinus();
if(checkConvergence(dof,energy,nModes)) {
fprintf(results,"Converged at iteration %d\n",iter);
fprintf(results,"Excited-State VSCF Energy is: %.8f\n", energy*219474.6313708);
excitedEnergies[z+1] = energy;
break;
} else {
prevEnergy = energy;
}
if(iter == 99) {
print(results,"VSCF failed to converge.\n");
excitedEnergies[z+1] = energy;
}}
}

///////////////////////////////////////DIPOLE CALCULATIONS/////////////////////////////



  //Print out all the transition frequencies
  print(results,"************************************************************************\n");
  print(results," VSCF ground-state energy: (cm^-1) \n");
  fprintf(results," % -15.4f \n", excitedEnergies[0]*(219474.6313708));
  print(results," \n");
  print(results," Transitions: (cm^-1) \n");
  print(results," \n");
  print(results,"  Harmonic        VSCF            Intensity(km/mol)\n");
  for(int i=0; i<nModes ; i++) {
  fprintf(results," % -15.4f % -15.4f \n", freq[i],(excitedEnergies[i+1]-excitedEnergies[0])*(219474.6313708));
  }

  //DEALLOCATE
    delete[] freq;
    delete[] mass;
    for(int i=0 ; i<nModes ; i++) {
      delete[] effV[i];
      delete[] slices[i];
      delete dof[i];
    }
    delete[] effV;
    delete[] slices;
    delete[] dof;
    delete[] excitedEnergies;
  }
  return 0;
}

//=====================================OTHER METHODS==========================================
void readin(double* freq, double* mass, int N) {
  //read in frequencies 
  std::ifstream in("freq.dat",std::ios::in);
  if(!in) {
    printf("Error: freq.dat could not be opened\n");
    exit(0);
  }

  for(int i=0 ; i<N ; i++) in >> freq[i];
  in.close();
  
  //read in reduced masses
  in.open("rmass.dat");
  if(!in) {
    printf("Error: rmass.dat could not be opened\n");
    exit(0);
  }
  for(int i=0 ; i<N ; i++) in >> mass[i];
  in.close();
} 


void readPot(std::string fileName, double* potential, int length, int potLength) {
  //read in potential and subtract out equilibrium energy
  std::ifstream in(fileName,std::ios::in);
  if(!in) {
    printf("Error: %s could not be opened\n",fileName.c_str());
    exit(0);
  }
  double value;
  int index = 0;
  while(in >> value) {
    if(index < potLength) {
      potential[index++] = value;
    } else {
      printf("%s is too long.\n",fileName.c_str());
      exit(0);
    } 
  }
  if(index < potLength) {
    printf("%s is too short.\n",fileName.c_str());
    exit(0);
  }
  in.close();

  double eqEnergy = potential[(length-1)/2];
  for(int i=0 ; i<potLength ; i++) {
    potential[i] -= eqEnergy;
  }
} 

bool checkConvergence(Mode** dof, double energy, int nModes) {
  double diff = 0.0;
  for(int i=0 ; i<nModes ; i++) {
    double temp = dof[i]->computeMaxDiff();
    if(temp > diff)
      diff = temp;
  }
  return (diff < 1.0E-5) && (fabs(energy-prevEnergy)*219474.6313708 <0.5);
}

int fact(int n) {
  if(n==0 || n==1)
    return 1;
  else
    return n*fact(n-1);
}

void print(FILE* script, std::string line) {
  fprintf(script,line.c_str());
}

int lengthCheck(int nMode, int coupling) {
  int num = nMode;
  int couplingFac = fact(coupling);
  for(int i=1; i<coupling ; i++)
    num *= (nMode-i);
  return num/couplingFac;
}
