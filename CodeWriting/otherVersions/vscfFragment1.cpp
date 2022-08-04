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
      readPot("Dx.dat",Dx,(int) pow(nPoints,couplingDegree),expectedLength);
      readPot("Dy.dat",Dy,(int) pow(nPoints,couplingDegree),expectedLength);
      readPot("Dz.dat",Dz,(int) pow(nPoints,couplingDegree),expectedLength);
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
