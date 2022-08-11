#ifndef AUX_H
#define AUX_H
#include <stdlib.h>
#include <stdio.h>
#include <mkl.h>

double* AllocDouble(int size);
void ABmult(double* C, double* A, double* B, int M, int N, int K, int LDC, int LDA, int LDB, int Type);
void diagonalize(double* pEigvecs, double* pEigvals, int n);
void linsolver(double* A, double* B, int N);

#endif
