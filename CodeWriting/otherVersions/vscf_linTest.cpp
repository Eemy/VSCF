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
  //SET-UP
    int nModes = atoi(argv[1]);
    int nPoints = atoi(argv[2]);
    int couplingDegree = atoi(argv[3]);
    int iter = 100;
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

    double* actualDx;
    double* actualDy;
    double* actualDz;
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
//      readPot("V.dat",V,(int) pow(nPoints,couplingDegree),expectedLength);

      //Read in the dipoles
//      Dx = new double[expectedLength];
//      Dy = new double[expectedLength];
//      Dz = new double[expectedLength];
//
      actualDx = new double[expectedLength];
      actualDy = new double[expectedLength];
      actualDz = new double[expectedLength];
      readPot("Dx.dat",actualDx,(int) pow(nPoints,couplingDegree),expectedLength);
      readPot("Dy.dat",actualDy,(int) pow(nPoints,couplingDegree),expectedLength);
      readPot("Dz.dat",actualDz,(int) pow(nPoints,couplingDegree),expectedLength);
      for(int i=0 ; i<expectedLength ; i++) {
        actualDx[i]*=debye_to_ea0;
        actualDy[i]*=debye_to_ea0;
        actualDz[i]*=debye_to_ea0;
      }     

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
    double* intensityComponents = new double[3*nModes];
    double* intensities = new double[nModes];
    FILE *results = fopen("eemVSCF.dat","w");
//////////////////
  //FILE *Vprint = fopen("V.datPrint","w");
  int pairIndex = 0;
  int testLength = (int) pow(nPoints,couplingDegree);
  for(int a=0 ; a<nModes ; a++) {
    for(int b=a+1 ; b<nModes ; b++) {
      for(int c=0 ; c<nPoints ; c++) {
        for(int d=0 ; d<nPoints ; d++) {
          double harmTemp = 0.5*dof[a]->getMass()*
                                    dof[a]->getOmega()*dof[a]->getOmega()*
                                    dof[a]->getPoint(c)*dof[a]->getPoint(c)+
                                    0.5*dof[b]->getMass()*
                                    dof[b]->getOmega()*dof[b]->getOmega()*
                                    dof[b]->getPoint(d)*dof[b]->getPoint(d);
          V[pairIndex*testLength+c*nPoints+d] = harmTemp;
          //fprintf(Vprint," % -15.8f \n", harmTemp);
        }
      }     
      pairIndex++;
    }
  }
//////////////////
//=====================================================================================
Potential pot(V,1,couplingDegree,nModes,nPoints,expectedLength,dof);
Potential dx(actualDx,1,couplingDegree,nModes,nPoints,expectedLength,dof);
Potential dy(actualDy,1,couplingDegree,nModes,nPoints,expectedLength,dof);
Potential dz(actualDz,1,couplingDegree,nModes,nPoints,expectedLength,dof);

double** slices = pot.get1DSlices();
double** dxslices = dx.get1DSlices();
double** dyslices = dy.get1DSlices();
double** dzslices = dz.get1DSlices();

/////////////////////////////////////////////////////////
//Test: Finite Difference Derivative of dipole
double *dipDerivX = new double[nModes];
double *dipDerivY = new double[nModes];
double *dipDerivZ = new double[nModes];

for(int i =0 ; i<nModes ; i++)  {
 dipDerivX[i] = (dxslices[i][6]-dxslices[i][4])/(dof[i]->getPoint(6)-dof[i]->getPoint(4)); 
 dipDerivY[i] = (dyslices[i][6]-dyslices[i][4])/(dof[i]->getPoint(6)-dof[i]->getPoint(4)); 
 dipDerivZ[i] = (dzslices[i][6]-dzslices[i][4])/(dof[i]->getPoint(6)-dof[i]->getPoint(4)); 
}

//allocate hessian matrix
//for each unique pair of modes
/*std::vector<double> mixedX;
std::vector<double> mixedY;
std::vector<double> mixedZ;
pairIndex = 0;
for(int i=0 ; i<nModes ; i++) {
  for(int j=i+1 ; j<nModes ; j++) {
    //DX
    double xpoint1 = actualDx[pairIndex*testLength+(testLength-1)/2+nPoints+1]; 
    double xpoint2 = actualDx[pairIndex*testLength+(testLength-1)/2+nPoints-1];
    double xpoint3 = actualDx[pairIndex*testLength+(testLength-1)/2-nPoints+1];
    double xpoint4 = actualDx[pairIndex*testLength+(testLength-1)/2-nPoints-1];

      //temp1: (+i):(j+1)-(j-1)/2dj
    double temp1 = (xpoint1-xpoint2)/(dof[j]->getPoint(6)-dof[j]->getPoint(4));
      //temp2: (-i):(j+1)-(j-1)/2dj
    double temp2 = (xpoint3-xpoint4)/(dof[j]->getPoint(6)-dof[j]->getPoint(4));
      //temp1-temp2/2di
    double doubleDeriv = (temp1-temp2)/(dof[i]->getPoint(6)-dof[i]->getPoint(4));
    mixedX.push_back(doubleDeriv);

    //DY
    double ypoint1 = actualDy[pairIndex*testLength+(testLength-1)/2+nPoints+1]; 
    double ypoint2 = actualDy[pairIndex*testLength+(testLength-1)/2+nPoints-1];
    double ypoint3 = actualDy[pairIndex*testLength+(testLength-1)/2-nPoints+1];
    double ypoint4 = actualDy[pairIndex*testLength+(testLength-1)/2-nPoints-1];

      //temp1: (+i):(j+1)-(j-1)/2dj
    temp1 = (ypoint1-ypoint2)/(dof[j]->getPoint(6)-dof[j]->getPoint(4));
      //temp2: (-i):(j+1)-(j-1)/2dj
    temp2 = (ypoint3-ypoint4)/(dof[j]->getPoint(6)-dof[j]->getPoint(4));
      //temp1-temp2/2di
    doubleDeriv = (temp1-temp2)/(dof[i]->getPoint(6)-dof[i]->getPoint(4));
    mixedY.push_back(doubleDeriv);

    //DZ
    double zpoint1 = actualDz[pairIndex*testLength+(testLength-1)/2+nPoints+1]; 
    double zpoint2 = actualDz[pairIndex*testLength+(testLength-1)/2+nPoints-1];
    double zpoint3 = actualDz[pairIndex*testLength+(testLength-1)/2-nPoints+1];
    double zpoint4 = actualDz[pairIndex*testLength+(testLength-1)/2-nPoints-1];

      //temp1: (+i):(j+1)-(j-1)/2dj
    temp1 = (zpoint1-zpoint2)/(dof[j]->getPoint(6)-dof[j]->getPoint(4));
      //temp2: (-i):(j+1)-(j-1)/2dj
    temp2 = (zpoint3-zpoint4)/(dof[j]->getPoint(6)-dof[j]->getPoint(4));
      //temp1-temp2/2di
    doubleDeriv = (temp1-temp2)/(dof[i]->getPoint(6)-dof[i]->getPoint(4));
    mixedZ.push_back(doubleDeriv);

    //printf("%.8f\n",mixedX[pairIndex]);
    //printf("%.8f\n",mixedY[pairIndex]);
    //printf("%.8f\n",mixedZ[pairIndex]);

    pairIndex++;
  }
}
*/
//Print derivs and points
//for(int i=0 ; i<nModes ; i++) {
//  printf("Mode %i DerivX: %.8f LeftPt, LeftVal: %.8f, %.8f RightPt, RightVal: %.8f, %.8f\n",i,dipDerivX[i],dof[i]->getPoint(4),dxslices[i][4],dof[i]->getPoint(6),dxslices[i][6]); 
//  printf("Mode %i DerivY: %.8f LeftPt, LeftVal: %.8f, %.8f RightPt, RightVal: %.8f, %.8f\n",i,dipDerivY[i],dof[i]->getPoint(4),dyslices[i][4],dof[i]->getPoint(6),dyslices[i][6]); 
//  printf("Mode %i DerivZ: %.8f LeftPt, LeftVal: %.8f, %.8f RightPt, RightVal: %.8f, %.8f\n",i,dipDerivZ[i],dof[i]->getPoint(4),dzslices[i][4],dof[i]->getPoint(6),dzslices[i][6]); 
//}

//Build Dipole Surfaces
  FILE *DxPrint = fopen("Dx.datPrint","w");
  FILE *DyPrint = fopen("Dy.datPrint","w");
  FILE *DzPrint = fopen("Dz.datPrint","w");

  Dx = new double[expectedLength];
  Dy = new double[expectedLength];
  Dz = new double[expectedLength];

  pairIndex = 0;
  testLength = (int) pow(nPoints,couplingDegree);
  for(int a=0 ; a<nModes ; a++) {
    for(int b=a+1 ; b<nModes ; b++) {
      for(int c=0 ; c<nPoints ; c++) {
        for(int d=0 ; d<nPoints ; d++) {
//          Dx[pairIndex*testLength+c*nPoints+d] = dipDerivX[a]+dipDerivX[b];
//          Dy[pairIndex*testLength+c*nPoints+d] = dipDerivY[a]+dipDerivY[b];
//          Dz[pairIndex*testLength+c*nPoints+d] = dipDerivZ[a]+dipDerivZ[b];

//          double dxTemp = dipDerivX[a]*dof[a]->getPoint(c)+dipDerivX[b]*dof[b]->getPoint(d)+mixedX[pairIndex]*dof[a]->getPoint(c)*dof[b]->getPoint(d);
//          double dyTemp = dipDerivY[a]*dof[a]->getPoint(c)+dipDerivY[b]*dof[b]->getPoint(d)+mixedY[pairIndex]*dof[a]->getPoint(c)*dof[b]->getPoint(d);
//          double dzTemp = dipDerivZ[a]*dof[a]->getPoint(c)+dipDerivZ[b]*dof[b]->getPoint(d)+mixedZ[pairIndex]*dof[a]->getPoint(c)*dof[b]->getPoint(d);

          double dxTemp = dipDerivX[a]*dof[a]->getPoint(c)+dipDerivX[b]*dof[b]->getPoint(d)+dof[a]->getPoint(c)*dof[b]->getPoint(d);//*dof[b]->getPoint(d);
          double dyTemp = dipDerivY[a]*dof[a]->getPoint(c)+dipDerivY[b]*dof[b]->getPoint(d)+dof[a]->getPoint(c)*dof[b]->getPoint(d);//*dof[b]->getPoint(d);
          double dzTemp = dipDerivZ[a]*dof[a]->getPoint(c)+dipDerivZ[b]*dof[b]->getPoint(d)+dof[a]->getPoint(c)*dof[b]->getPoint(d);//*dof[b]->getPoint(d);

          Dx[pairIndex*testLength+c*nPoints+d] = dxTemp;
          Dy[pairIndex*testLength+c*nPoints+d] = dyTemp;    
          Dz[pairIndex*testLength+c*nPoints+d] = dzTemp;

          fprintf(DxPrint," % -15.8f \n", dxTemp);
          fprintf(DyPrint," % -15.8f \n", dyTemp);
          fprintf(DzPrint," % -15.8f \n", dzTemp);
        }
      }     
      pairIndex++;
    }
  }

Potential dxLin(Dx,1,couplingDegree,nModes,nPoints,expectedLength,dof);
Potential dyLin(Dy,1,couplingDegree,nModes,nPoints,expectedLength,dof);
Potential dzLin(Dz,1,couplingDegree,nModes,nPoints,expectedLength,dof);

double** dxLinslices = dxLin.get1DSlices();
double** dyLinslices = dyLin.get1DSlices();
double** dzLinslices = dzLin.get1DSlices();
////////////////////////////////////////////////////////

//Prepare: eigensolver on pure 1D slices for each mode
for(int i = 0 ; i< nModes ; i++) {
prevEnergy += solver.solveMode(dof[i],slices[i],0);
}
//Save ground state
for(int i=0 ; i<nModes ; i++) 
 dof[i]->setGroundState(); 

for(int iter = 1 ; iter< 100 ; iter++) {
int counter = 0;
for(int i = 0 ; i< nModes ; i++) {
for(int j = 0 ; j< nPoints ; j++) {
effV[i][j] = slices[i][j];
}
}
for(int a = 0 ; a< nModes ; a++) {
for(int b = a+1 ; b< nModes ; b++) {
for(int c = 0 ; c< nPoints ; c++) {
effV[a][c] += pot.integralDriver(counter, 0, c);
effV[b][c] += pot.integralDriver(counter, 1, c);
}
counter++;
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
//Save ground state
for(int i=0 ; i<nModes ; i++) 
 dof[i]->setGroundState(); 
}
if(iter == 99) {
print(results,"VSCF failed to converge.\n");
excitedEnergies[0] = energy;
}}
/////////////////End Ground-State VSCF/////////////////////


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
for(int c = 0 ; c< nPoints ; c++) {
effV[a][c] += pot.integralDriver(counter, 0, c);
effV[b][c] += pot.integralDriver(counter, 1, c);
}
counter++;
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
dof[z]->setExcitedState();
}


///////////////////////////////////////DIPOLE CALCULATIONS/////////////////////////////
//1D Integrals
for(int i=0 ; i<nModes ; i++) {
 // double factor = 1/(sqrt(2*dof[i]->getMass()*dof[i]->getOmega()));
 // double overlap = dof[i]->getOverlapEG();
//  printf("Factor: %.8f\n",factor);
//  printf("1Dx: %.8f, actual: %.8f\n",dx.integrateSlice(dof[i],dxLinslices[i],true),dipDerivX[i]*factor);
//  printf("1Dy: %.8f, actual: %.8f\n",dy.integrateSlice(dof[i],dyLinslices[i],true),dipDerivY[i]*factor);
//  printf("1Dz: %.8f, actual: %.8f\n",dz.integrateSlice(dof[i],dzLinslices[i],true),dipDerivZ[i]*factor);

  intensityComponents[3*i] += dx.integrateSlice(dof[i],dxLinslices[i],true);
  intensityComponents[3*i+1] += dy.integrateSlice(dof[i],dyLinslices[i],true);
  intensityComponents[3*i+2] += dz.integrateSlice(dof[i],dzLinslices[i],true);
}

//Integrals beyond 1D
int counter = 0;
for(int a = 0 ; a< nModes ; a++) {
  for(int b = a+1 ; b< nModes ; b++) {
    //Intensity Contribution for Mode A 
    //<0|D(q1)|0> <0|1>
    double overlapA = dof[a]->getOverlapEG();
    intensityComponents[3*a] += dx.integrateSlice(dof[b],dxLinslices[b],false)*overlapA;
    intensityComponents[3*a+1] += dy.integrateSlice(dof[b],dyLinslices[b],false)*overlapA;
    intensityComponents[3*a+2] += dz.integrateSlice(dof[b],dzLinslices[b],false)*overlapA;

    //<10|~D(q1,q2)|00>
    printf("2DX Pair %i,%i, ModeA: %.8f, Expected: %.8f\n",a,b,dxLin.getDipole(counter,0),
  sqrt(1/(2*dof[a]->getMass()*dof[a]->getOmega()))*1/(2*dof[b]->getMass()*dof[b]->getOmega()));
    printf("2DY, ModeA: %.8f\n",dyLin.getDipole(counter,0));
    printf("2DZ, ModeA: %.8f\n",dzLin.getDipole(counter,0));
    intensityComponents[3*a] += dxLin.getDipole(counter,0);
    intensityComponents[3*a+1] += dyLin.getDipole(counter,0);
    intensityComponents[3*a+2] += dzLin.getDipole(counter,0);

    //Intensity Contribution for Mode B
    double overlapB = dof[b]->getOverlapEG();
    intensityComponents[3*b] += dx.integrateSlice(dof[a],dxLinslices[a],false)*overlapB;
    intensityComponents[3*b+1] += dy.integrateSlice(dof[a],dyLinslices[a],false)*overlapB;
    intensityComponents[3*b+2] += dz.integrateSlice(dof[a],dzLinslices[a],false)*overlapB;

    
    printf("2DX, ModeB: %.8f\n",dxLin.getDipole(counter,1)); //should be 0
    printf("2DY, ModeB: %.8f\n",dyLin.getDipole(counter,1));
    printf("2DZ, ModeB: %.8f\n",dzLin.getDipole(counter,1));
    intensityComponents[3*b] += dxLin.getDipole(counter,1);
    intensityComponents[3*b+1] += dyLin.getDipole(counter,1);
    intensityComponents[3*b+2] += dzLin.getDipole(counter,1);
    counter++;
  }
}
for(int i=0 ; i<nModes ; i++) {
  for(int j=0 ; j<3 ; j++) {
//    printf("Component %i: %.8f\n",j,intensityComponents[3*i+j]);
    intensities[i] += intensityComponents[3*i+j]*intensityComponents[3*i+j]; 
  }
//  printf("Before Prefactor: %.8f\n",intensities[i]);
//  printf("Freq: %.8f\n",excitedEnergies[i+1]-excitedEnergies[0]);
  intensities[i] *= (excitedEnergies[i+1]-excitedEnergies[0])*2.0*pi*Na/(3.0*c_in_au*c_in_au)*a0_to_cm*1.0E-5; //to km/mol
}
//printf("Prefactor: %.8f\n",2.0*pi*Na/(3.0*c_in_au*c_in_au)*a0_to_cm*1.0E-5);

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
//  fprintf(results," % -15.4f % -15.4f \n", freq[i],(excitedEnergies[i+1]-excitedEnergies[0])*(219474.6313708));
  }

  //DEALLOCATE
    delete[] freq;
    delete[] mass;
    for(int i=0 ; i<nModes ; i++) {
      delete[] effV[i];
      delete[] slices[i];
      delete[] dxslices[i];
      delete[] dyslices[i];
      delete[] dzslices[i];
      delete dof[i];
    }
    delete[] effV;
    delete[] slices;
    delete[] dxslices;
    delete[] dyslices;
    delete[] dzslices;
    delete[] dof;
    delete[] excitedEnergies;
    delete[] intensities;
    delete[] intensityComponents;
/////////////test dip////////////
//delete[] groundHarm;
//delete[] excitedHarm;
/////////////////////////////////
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
