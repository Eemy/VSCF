#ifndef MP2CORR_H
#define MP2CORR_H
#include <vector>
#include <utility>
#include <stdbool.h>
using std::pair;
using std::vector;

class Mode;
class Potential;
class Mp2Corr {
  public:
    Mp2Corr(vector<Mode*>&, vector<Potential*>&, vector<vector<vector<int>>>&); 
    ~Mp2Corr();
    void getSecondOrderCorr(vector<double>&);
    void getSecondOrderCorr(vector<double>&,double*);
    void calculateIntegrals(vector<pair<vector<int>,vector<int>>>,vector<double>);
    void clear();
  private:
    int maxState;
    int minState;
    int numStates;
    int nModes;
    int maxDim;
    int numPsi;

    vector<Mode*> dof;
    vector<Potential*> pot;
    vector<vector<vector<int>>> tuples;

    vector<double> singles;
    vector<double> singlesDenom;
    vector<vector<double>> integrals;
    vector<vector<double>> denominators;

    void fillCorrectionMatrices(int,int,int,int,int,vector<int>,vector<double>&);
    int getIndex(vector<int>,int,int);
    bool integralIsNonZero(vector<int>,vector<int>); 
    bool brillouin();
};

#endif
