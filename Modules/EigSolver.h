#ifndef EIGSOLVER_H
#define EIGSOLVER_H
#include <vector>

class Mode;
class EigSolver {
    int nPoints;
    int nBasis;
    double* T;
    void buildTmat();
  public:
    EigSolver(int);
    ~EigSolver();
    double solveMode(Mode*, std::vector<double>, int);
};

#endif
