#ifndef MODE_H
#define MODE_H

class Mode {
  public:
    Mode(double, double, int);
    ~Mode();
    void setGroundState(double*); 
    void updateWaveFcn(double*);
    double* getWaveFcn();
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
    double m;
    int nPoints;
    int nBasis;
    double* waveFcn;
    double* oldWaveFcn;
    double* points;
    double* weights;
    double* hermiteEval;
    double* norm;
    
    double* groundState;
};

#endif
