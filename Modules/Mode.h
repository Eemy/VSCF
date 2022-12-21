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
    void useVSCFStates(bool);
//    void updateWaveFcn(double*);
    void updateAllPsi_AllE(double*,double*);
    double* getModal(int);
    double getEModal(int); 
    double getEModal();
//    double* getWaveFcn();
//    double* getGState(); 
//    double* getEState(); 
    double getOverlapEG();
    double computeMaxDiff();
//    double getAlpha();
    double getOmega();    
//    double getMass();
    int getNPoints();
    double getWeight(int);
    double getHerm(int, int);
    double getNorm(int);
    void setHarmonic();
    void setAnharmonic();
    double getIntegralComponent(int);
//    double getPoint(int);
    double *getDensity();
    void setBra(int);
    void setKet(int);
    int getBra();
    int getKet();
    void saveErrorVec(double*,int);
    int getNumErrorVecs();
    double dotErrorVecs(int, int);
    void resetSubspace();
    void extrapolateDensity(double *);
//    void diis(double*,double*,int);
    double getDIISError();
 private:
    void setMaxElement(double*);
    void updateDensity();

    double omega;
    int nPoints;
    int nBasis;
//    double* waveFcn;
//    double* oldWaveFcn;
    double* energies;
    double* waveAll;
    double* waveAll_prev;
    double* tempAll;

    double* points;
    double* weights;
    double* hermiteEval;
    double* norm;
    
    double* vscfPsi;
    bool vscfStates; 
    bool harm;
    int bra;
    int ket;
  
    double maxDiisError;
    int conv;
//    std::vector<double*> Fsave;
    std::vector<double*> Dsave;
    std::vector<double*> Esave; 
    double* density;
    int diis_subspace;
};

#endif
