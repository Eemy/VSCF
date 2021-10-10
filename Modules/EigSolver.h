#ifndef EIGSOLVER_H
#define EIGSOLVER_H

class Mode;
class EigSolver {
    int nBasis;
    double* T;
    void buildTmat();
  public:
    EigSolver(int);
    ~EigSolver();
    double solveMode(Mode*, double*, int);
};

#endif
