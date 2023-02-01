#include <cstdio>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include "../UtilityFunc/ghQuad.h"

#define pi  3.1415926535897932384626433832795

int main(int argc, char* argv[]) {
  if(argc < 4) { 
    printf("Error: <nPoints> <coefficient File> <1=gh, 2=even-space>\n"); 
    exit(0);
  }

  int nPoints = atoi(argv[1]);
  std::string file = argv[2];
  int gridChoice = atoi(argv[3]);

  //Read basis fcn coeff file
  std::ifstream in(file,std::ios::in);  
  if(!in) {
    printf("Error: %s could not be opened.\n",file.c_str());
    exit(0);
  }
  std::vector<double> coeff; 
  double val = 0.0;
  while(in >> val) {
    coeff.push_back(val);
  } 
  int nBasis = coeff.size(); 

  //Generate HO wavefunctions
  double *points = new double[nPoints];
  double *hermiteEval = new double[nPoints*nBasis];
  double *norm = new double[nBasis];
  //for ghGrid:
  if(gridChoice == 1) {
    double *weights = new double[nPoints];   
    gauher(points,weights,nPoints);
    delete[] weights;
  //for evenly-spaced Grid:
  } else if(gridChoice == 2) {
    double min = -10.0;
    double max = 10.0;    
    double delta = (max-min)/(nPoints-1);
    for(int i=0 ; i<nPoints ; i++) points[i] = min+i*delta;
  }

  for(int i=0 ; i<nBasis ; i++) {
    for(int j=0 ; j<nPoints ; j++) {
      hermiteEval[i*nPoints+j] = hermite(i,points[j]);
    }
    norm[i] = 1/(sqrt(pow(2.0,i)*factorial(i)))*pow(1/pi,0.25);
  }

  //Calculate value at each gridpoint
  std::vector<double> vals(nPoints);
  for(int i=0 ; i<nPoints; i++) {
    for(int j=0 ; j<nBasis ; j++) {
      vals[i] += coeff[j]*norm[j]*hermiteEval[j*nPoints+i]
                *exp(-(points[i]*points[i]/2));
    }
    printf("%.8f\t%.8f\n",points[i],vals[i]);
  }

  delete[] norm;
  delete[] points;
  delete[] hermiteEval; 
}
