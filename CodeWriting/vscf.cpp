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

void readin(std::vector<Mode*>& dof, int N, int nPoints);
void readPot(std::string, double* potential, int length, int potLength);
bool checkConvergence(Mode** dof, double energy, int nModes);
int fact(int n);
void print(FILE* script, std::string line); 
int lengthCheck(int nMode, int coupling);

double prevEnergy = 0.0;

int main(int argc, char* argv[]) {
  const int defaultLength = 8;
  int maxIter = 100;

  if(argc < defaultLength) { 
    printf("Error: <Nmodes> <Nquad> <EnergyFile> <DipoleXFile> <DipoleYFile> <DipoleZFile> <CouplingDegree> [<EnergyFile> <Dx> <Dy> <Dz> <dim> ...]\n");
    exit(0);
  }

  //SET-UP AND READ IN ARGS
  int nModes = atoi(argv[1]); //arg1
  int nPoints = atoi(argv[2]); //arg2
  std::vector<std::string> potFileNames;
  potFileNames.push_back(argv[3]); //arg3
  std::vector<std::string> dipFileNames;
  dipFileNames.push_back(argv[4]); //arg4
  dipFileNames.push_back(argv[5]); //arg5
  dipFileNames.push_back(argv[6]); //arg6
  std::vector<int> potDims;
  potDims.push_back(atoi(argv[7])); //arg7

/////////////////////////Create Mode, EigSolver, Potential Objects////////////////////////
  std::vector<Mode*> dof;
  readin(dof,nModes,nPoints); 
  EigSolver solver(nPoints);

  std::vector<Potential*> pot;
  std::vector<Potential*> dip;
  //For first set of V and D
  pot.push_back(new Potential(potFileNames[0],1,potDims[0],nModes,nPoints,dof));
  for(int i=0 ; i<3 ; i++) 
    dip.push_back(new Potential(dipFileNames[i],1,potDims[0],nModes,nPoints,dof));

  //For subsequent sets of V and D (5 represents number of args one set occupies in argv)
  if ((argc-defaultLength)%5 == 0 && (argc-defaultLength) > 0) {
    int len = (argc-defaultLength)/5;
    for(int i=1 ; i<=len ; i++) {
      potFileNames.push_back(argv[defaultLength+5*(i-1)]);
      potDims.push_back(atoi(argv[defaultLength+5*(i-1)+4]));
      pot.push_back(new Potential(potFileNames[i],potDims[i-1]+1,potDims[i],nModes,nPoints,dof)); 
      for(int j=1 ; j<=3 ; j++) {
        dipFileNames.push_back(argv[defaultLength+5*(i-1)+j]);
        dip.push_back(new Potential(dipFileNames[i*3+(j-1)],potDims[i-1]+1,potDims[i],nModes,nPoints,dof));
      }
    } 
  } else {
    printf("The number of args is invalid. Check your input and try again.\n");
    exit(0);
  }

  double** slices = pot[0]->get1DSlices();
  std::vector<double**> dipSlices;
  dipSlices.push_back(dip[0]->get1DSlices());
  dipSlices.push_back(dip[1]->get1DSlices());
  dipSlices.push_back(dip[2]->get1DSlices());

  //Allocate a few additional arrays for the calculations
  double** effV = new double*[nModes];
  for(int i=0 ; i<nModes ; i++) {
    effV[i] = new double[nPoints];
  }
  double* excitedEnergies = new double[nModes+1];
  double* intensityComponents = new double[3*nModes];
  double* intensities = new double[nModes];
  std::vector<double> overlaps;

  //Open results file once all set-up is completed
  FILE *results = fopen("eemVSCF.dat","w");

/*
  readPot("V.dat",V,(int) pow(nPoints,couplingDegree),expectedLength);
  readPot("D3x.dat",Dx,(int) pow(nPoints,couplingDegree),expectedLength);
  readPot("D3y.dat",Dy,(int) pow(nPoints,couplingDegree),expectedLength);
  readPot("D3z.dat",Dz,(int) pow(nPoints,couplingDegree),expectedLength);
  for(int i=0 ; i<expectedLength ; i++) {
    Dx[i]*=debye_to_ea0;
    Dy[i]*=debye_to_ea0;
    Dz[i]*=debye_to_ea0;
 }
*/
//====================================Begin VSCF============================================
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
  delete[] excitedEnergies;
  delete[] intensities;
  delete[] intensityComponents;
  for(int i=0 ; i<dip.size() ; i++) delete dip[i];
  for(int i=0 ; i<pot.size() ; i++) delete pot[i];

  return 0;
}

//=====================================OTHER METHODS==========================================
void readin(std::vector<Mode*>& dof, int N, int nPoints) {
  std::vector<double> freq;
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
    dof.push_back(new Mode(freq[i],nPoints));
  }    
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
