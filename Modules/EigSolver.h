#ifndef EIGSOLVER_H
#define EIGSOLVER_H
#include <vector>

class Mode;
class EigSolver {
    int nPoints;
    int nBasis;
    int conv;
    double* T;
    void buildTmat();
  public:
    EigSolver(int,int);
    ~EigSolver();
    double solveMode(Mode*, std::vector<double>, int, int);
//    double solveMode(Mode*, std::vector<double>, int, int, int);
//    void diis(std::vector<Mode*>);
    void diis(std::vector<Mode*>,int);
//    void diis(Mode*);
};

#endif
