#ifndef POTENTIAL_H
#define POTENTIAL_H
#include <vector>
#include <string>
#include <unordered_map>
#include <stdbool.h>

class cpot;
class Mode;
class Potential {
  public:
    Potential(std::string, int, int, int, int);
    ~Potential();
    std::vector<std::vector<double>> get1DSlices();
    double integralDriver(int, int, int);
    double integrateTuple(int, bool);
    double getDipole(int, int);
    double integrateSlice(Mode*, int);
    std::vector<std::vector<int>> readPot(std::vector<Mode*>&, bool);
    int nPoints;
    int nModes;
  private:
    class Tuplet {
      private:
        class SubTup {
          public: 
            SubTup(std::vector<int>, std::vector<Mode*>, std::vector<double>);
            ~SubTup();
            std::vector<int> subspaceSet;
            std::vector<Mode*> modeSet;
            std::vector<double> potSet;
        };
      public:
        Tuplet(int, int, int, std::vector<double>, std::vector<Mode*>, std::vector<int>);
        ~Tuplet();
        std::vector<Mode*> modeSubset;
        std::vector<int> indices;
        std::vector<std::vector<SubTup*>> subTuplets;
        std::vector<double> pot;
        int dim;
        int minSubDim;
        int nPoints;
        int potLength;
        std::vector<int> subspaces;
        double getEffVIntegral(int, int);
        double getTotalIntegral(bool);
        cpot *couplingRemover;
        void setUpDriver(std::unordered_map<std::string,int>&);

        void getSlice(int, int, int, std::vector<double>&, std::vector<double>&, std::vector<int>&);
        double getIntegral(int, int, int, int, std::vector<Mode*>&, std::vector<double>&, int*);
        void fillPotential(int, int, std::vector<double>&, std::vector<double>&, int, int*, std::vector<int>&);
        void setUp(int, std::vector<double>&, std::unordered_map<std::string,int>&, std::vector<int>&, std::vector<int>&, std::vector<Mode*>&);
    };
  private: 
    std::string file;
    int minDim;
    int dim;
    int potLength;
    std::vector<std::vector<double>> slices;
    std::vector<int> subspaces;
    std::vector<Tuplet*> tuplets;
    //std::unordered_map<std::string,Tuplet*> tuplets;
    std::unordered_map<std::string,int> tupletChecker;
};

#endif
