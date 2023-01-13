#ifndef MP2CORR_H
#define MP2CORR_H
#include <vector>

class Mp2Corr {
  public:
  private:
    int maxState;
    int minState;
    int numStates;
    int nModes;

    std::vector<Mode*> dof;
    std::vector<Potential*> pot;

    std::vector<double> singles;
    std::vector<double> singlesDenom;
    std::vector<std::vector<double>> integrals;
    std::vector<std::vector<double>> denominators;
};

#endif
