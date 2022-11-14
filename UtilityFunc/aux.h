#ifndef AUX_H
#define AUX_H
#include <stdlib.h>
#include <stdio.h>
#include <mkl.h>

double* AllocDouble(int size);
void ABmult(double* C, double* A, double* B, int M, int N, int K, int LDC, int LDA, int LDB, int Type);
void diagonalize(double* pEigvecs, double* pEigvals, int n);
void linsolver(double* A, double* B, int N);
void printmat(double* mat, int n, int m, double scale);
int tupleIndexDriver(std::vector<int>& targetIndices, int nModes);
int tupleIndex(int& counter, int dim, int recursionLevel, std::vector<int>& indices, std::vector<int>& targetIndices, int nModes);

#endif
