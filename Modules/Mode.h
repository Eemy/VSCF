#ifndef MODE_H
#define MODE_H
#include <stdbool.h>
#include <vector>

class Mode {
  public:
    Mode(double, int, int);
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
    //double getPoint(int);
    double getDIISError();
    double *getDensity();
    void diis(double*,double*,int);
  private:
    void setMaxElement(double*);
    void updateDensity();

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

    int conv;
    std::vector<double*> Fsave;
    std::vector<double*> Esave; 
    double* density;
    double maxDiisError;
    int diis_subspace;
};

#endif
