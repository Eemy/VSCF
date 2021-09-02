#ifndef CPOT_H
#define CPOT_H
#include<string.h>
#include<iostream>
#include<vector>
using namespace std;

class cpot{ 

public:
//cpot(double*,int,int,int);
cpot(int);
~cpot();
//double get_coupling(double*,int,int,int);
void get_coupling(vector<double>&,int);
//double get2d(double*,int,int);
void get2d(vector<double>&);
//double get3d(double*,int,int);
void get3d(vector<double>&);
void get4d(vector<double>&);
void get5d(vector<double>&);
void get6d(vector<double>&);
int fact(int);
private:
//int degree;
int grid_size;
//int dof;

};
#endif

