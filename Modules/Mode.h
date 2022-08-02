#ifndef MODE_H
#define MODE_H
#include <stdbool.h>

class Mode {
  public:
    Mode(double, int);
    ~Mode();
    void setGroundState(); 
    void setExcitedState();
    void setExcited(bool);
    void updateWaveFcn(double*);
    double* getWaveFcn();
    double* getGState(); 
    double* getEState(); 
    double getOverlapEG();
    double computeMaxDiff();
    double getAlpha();
    double getOmega();    
    double getMass();
    int getNPoints();
    double getWeight(int);
    double getHerm(int, int);
    double getNorm(int);
    double getIntegralComponent(int);
    double getPoint(int);
  private:
    double alpha;
    double omega;
    int nPoints;
    int nBasis;
    double* waveFcn;
    double* oldWaveFcn;
    double* points;
    double* weights;
    double* hermiteEval;
    double* norm;
    
    double* groundState;
    double* excitedState;
    bool excited;
};

#endif
