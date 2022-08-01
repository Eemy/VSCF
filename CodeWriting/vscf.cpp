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

void readin(double* freq, double* mass, int N);
void readPot(std::string, double* potential, int length, int potLength);
bool checkConvergence(Mode** dof, double energy, int nModes);
int fact(int n);
void print(FILE* script, std::string line); 
int lengthCheck(int nMode, int coupling);

double prevEnergy = 0.0;

int main(int argc, char* argv[]) {
  if(argc!=5) printf("Error: <Nmodes> <Nquad> <couplingDegree> <#ofMethodsUsed>\n");
  else {
  //SET-UP AND READ IN ARGS
    int nModes = atoi(argv[1]);
    int nPoints = atoi(argv[2]);
    int couplingDegree = atoi(argv[3]);
    int maxIter = 100;
    int numMethods = atoi(argv[4]);    

    std::vector<int> potDims;
    std::vector<int> expectedLengths;
    int expectedLength = 0;
    std::vector<std::string> potFileNames;
    std::vector<double*> potentials;
    double *V;
    double* Dx;
    double* Dy;
    double* Dz;
    std::vector<double> overlaps;
//===================================Read in Potential(s) and Dipoles======================================
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

      //Read in the dipoles
      Dx = new double[expectedLength];
      Dy = new double[expectedLength];
      Dz = new double[expectedLength];
      readPot("D3x.dat",Dx,(int) pow(nPoints,couplingDegree),expectedLength);
      readPot("D3y.dat",Dy,(int) pow(nPoints,couplingDegree),expectedLength);
      readPot("D3z.dat",Dz,(int) pow(nPoints,couplingDegree),expectedLength);
      for(int i=0 ; i<expectedLength ; i++) {
        Dx[i]*=debye_to_ea0;
        Dy[i]*=debye_to_ea0;
        Dz[i]*=debye_to_ea0;
      }     
    }
//===================================A few more arrays...==========================================

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
    double* intensityComponents = new double[3*nModes];
    double* intensities = new double[nModes];
    //Open results file once all set-up is completed
    FILE *results = fopen("eemVSCF.dat","w");
//=================================================================================================
Potential pot(V,1,couplingDegree,nModes,nPoints,expectedLength,dof);
Potential dx(Dx,1,couplingDegree,nModes,nPoints,expectedLength,dof);
Potential dy(Dy,1,couplingDegree,nModes,nPoints,expectedLength,dof);
Potential dz(Dz,1,couplingDegree,nModes,nPoints,expectedLength,dof);

double** slices = pot.get1DSlices();

std::vector<double**> dipSlices;
dipSlices.push_back(dx.get1DSlices());
dipSlices.push_back(dy.get1DSlices());
dipSlices.push_back(dz.get1DSlices());
std::vector<Potential*> dipoles;
dipoles.push_back(&dx);
dipoles.push_back(&dy);
dipoles.push_back(&dz);

//Prepare: eigensolver on pure 1D slices for each mode
for(int i = 0 ; i< nModes ; i++) {
prevEnergy += solver.solveMode(dof[i],slices[i],0);
}
for(int iter = 1 ; iter< maxIter ; iter++) {
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
if(iter == maxIter-1) {
print(results,"Ground-State VSCF failed to converge.\n");
excitedEnergies[0] = energy;
}
}
///////End Ground-State VSCF/////////
for(int i = 0 ; i< nModes ; i++) {
dof[i]->setGroundState();
}
/////////Excited-State VSCF//////////
for(int z = 0 ; z< nModes ; z++) {
prevEnergy = 0.0;
for(int i = 0 ; i< nModes ; i++) {
if(i==z) {
prevEnergy += solver.solveMode(dof[i],slices[i],1);
} else {
 prevEnergy += solver.solveMode(dof[i],slices[i],0);
}
}
for(int iter = 1 ; iter< maxIter ; iter++) {
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
}
////////End Excited-State VSCF///////

for(int i = 0 ; i< nModes ; i++) {
dof[i]->updateWaveFcn(dof[i]->getGState());
overlaps.push_back(dof[i]->getOverlapEG());
}

//DIPOLE CALCULATIONS
for(int comp = 0 ; comp< 3 ; comp++) {
int counter = 0;
for(int a = 0 ; a< nModes ; a++) {
intensityComponents[3*a+comp] += dipoles[comp]->integrateSlice(dof[a],dipSlices[comp][a],true);
for(int b = a+1 ; b< nModes ; b++) {
intensityComponents[3*a+comp] += dipoles[comp]->integrateSlice(dof[b],dipSlices[comp][b],false)*overlaps[a];
intensityComponents[3*b+comp] += dipoles[comp]->integrateSlice(dof[a],dipSlices[comp][a],false)*overlaps[b];
for(int c = b+1 ; c< nModes ; c++) {
intensityComponents[3*a+comp] += dipoles[comp]->getDipole(counter,0);
intensityComponents[3*b+comp] += dipoles[comp]->getDipole(counter,1);
intensityComponents[3*c+comp] += dipoles[comp]->getDipole(counter,2);
counter++;
}
}
}
}
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
    delete[] freq;
    delete[] mass;
    for(int i=0 ; i<nModes ; i++) {
      delete[] effV[i];
      delete[] slices[i];
      delete[] dipSlices[0][i];
      delete[] dipSlices[1][i];
      delete[] dipSlices[2][i];
      delete dof[i];
    }
    delete[] effV;
    delete[] slices;
    delete[] dipSlices[0];
    delete[] dipSlices[1];
    delete[] dipSlices[2];
    delete[] dof;
    delete[] excitedEnergies;
    delete[] intensities;
    delete[] intensityComponents;
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
