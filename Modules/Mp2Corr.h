#ifndef MP2CORR_H
#define MP2CORR_H
#include <vector>
#include <pair>
using std::pair;
using std::vector;

class Mp2Corr {
  public:
    Mp2Corr(vector<Mode*>&, vector<Potential*>&, vector<vector<vector<int>>>&); 
    ~Mp2Corr();
    calculateIntegrals(vector<pair<vector<int>,vector<int>>>,vector<double>);
    clear();
  private:
    int maxState;
    int minState;
    int numStates;
    int nModes;
    int maxDim;

    vector<Mode*> dof;
    vector<Potential*> pot;
    vector<vector<vector<int>>> tuples;

    vector<double> singles;
    vector<double> singlesDenom;
    vector<vector<double>> integrals;
    vector<vector<double>> denominators;

    fillCorrectionMatrices(int,int,int,int,vector<int>,vector<double>&);
};

#endif
