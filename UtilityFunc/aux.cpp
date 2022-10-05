#include <cmath>
#include "aux.h"
//#include <cstdio>

using namespace std;


//=============================================
double* AllocDouble(int size)
{
  double* thearray = (double*) malloc(sizeof(double)*size);
  if(thearray==NULL)
  {
    printf("Allocation bombed.\n");
    printf(" size = %3.2f bits (%3.2f MB)\n",(double)size*sizeof(double), sizeof(double)*(double)size / (1024.0*1024.0));
    exit(1);
  }
  int i;
  for (i = 0; i < size; i++){
    thearray[i] = 0.0;
  }
  return thearray;
}
//=============================================
// AtimsB from Q-Chem. Fortran to C port

// https://www.ibm.com/docs/en/essl/6.3?topic=mos-sgemm-dgemm-cgemm-zgemm-combined-matrix-multiplication-addition-general-matrices-their-transposes-conjugate-transposes
//
  // C = alpha * op(A) * op(B) + beta * C
  //
  // Matrix A must always be stored in untransposed form.
  //
  // (transa, transb, l, n, m, alpha, a, lda, b, ldb, beta, c, ldc); 
  //////////////
  // transa == transb == 'N', 'T' or 'C'. The latter is complex conjugate a.k.a. Hermitian of a matrix.
  // l = int number of rows in matrix C
  // n = int number of columns in matrix C
  // m = number of columns in A and rows in B == dimension of "contraction".
  // alpha = scalar (should be taken care of via wrapper)
  // 
  // A == op(A). If 'N' -> (l*m) and LDA >= m. If 'T' || 'C' -> (m*l) and LDA >= l.
  // LDA = Leading dimension of A
  // B == op(B). If 'N' -> (m*n) and LDB >= n. If 'T' || 'C' -> (n*m) and LDB >= m.
  // LDB = leading dimension of B
  // beta = scalar (should be taken care of via wrapper)
  // 
  // C = is a (l*n) resulting matrix. 
  // LDC = leading dimension of C, must be at least int LDC >= n.
  // 
void ABmult(double* C, double* A, double* B, int M, int N, int K, int LDC, int LDA, int LDB, int Type){
  
  double alpha;
  double beta;
  alpha = 0.0;
  beta = 0.0;
  
  int jobtype;
  
  char transA;
  char transB;

  if (Type > 0){
    alpha = 1.0;
  }
  else{
    alpha = -1.0;
  }
  
  jobtype = fabs(Type);
  if (jobtype > 10){
    beta = 1.0;
  }
  else{
    beta = 0.0;
  }
  
  //printf("\n assigning transA and transB \n"); 
  transA = 'N'; // normal, not transposed matrix form
  transB = 'N'; // ditto
  //printf("\n assigned transA and transB \n"); 
  

  //printf("\n starting transA and B transpose check \n"); 
  jobtype %= 10;
  if (jobtype == 2 || jobtype == 4) transA = 'T';
  if (jobtype == 3 || jobtype == 4) transB = 'T';
  //printf("\n finished the check \n"); 
  
  //printf("\n starting the multiplier \n"); 
  dgemm(&transA, &transB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
  //printf("\n finished the multiplier \n"); 
  
  if (M < 0 || N < 0 || K < 0 || LDA < 0 || LDB < 0){
    printf("\n Invalid dimension in ABmult \n");
    exit(1);
  }
}

//=============================================
// matrix diagonalizer
void diagonalize(double* pEigvecs, double* pEigvals, int n){
  
  int i, j;
  int lda = n;
  int info;
  int lwork;
  double wkopt;
  double* work;

  // query and allocate the optimal working space
  lwork = -1; // initiates the workspace query
  dsyev("Vectors", "Upper", &n, pEigvecs, &lda, pEigvals, &wkopt, &lwork, &info);
  lwork = (int)wkopt;
  work = (double*)malloc(lwork*sizeof(double));
  
  // Solve eigenproblem
  dsyev("Vectors", "Upper", &n, pEigvecs, &lda, pEigvals, work, &lwork, &info);
  
  // check for crashes/convergence
  if (info > 0){
    printf("\n dsyev is GG \n");
    exit(1);
  }

  free((void*)work);
}
//=============================================
void linsolver(double* A, double* B, int N)
{

        //returns solution in B.
        int i;
        int NRS = 1;
        int LDA = N;
        int LDB = N;
        double Berr[1];
        double X[N];
        int IPIV[N];
        int IWork[N];
        double Ferr[1];
        double Work[4*N];
        double AF[LDA*N];
        double RCond, R, C;
        char fact = 'N';
        char trans= 'T';
        char equed = 'X';
        int info;

        dgesvx_(&fact,&trans,&N,&NRS,A,&LDA,AF,&LDA,IPIV,&equed,&R,&C,B,&LDB,X,&LDB,&RCond,Ferr,Berr,Work,IWork,&info);
        //printf("Rcond = % -10.16e \n",RCond); 
        for(i=0; i<N; i++) B[i] = X[i];
}

void printmat(double* mat, int n, int m)
{
int i,j;
for(i=0; i<n; i++){
  for(j=0; j<m; j++){
    printf("% -10.6e  ",mat[i*m+j]);
  }
 printf("\n");
 }
printf("\n");
}
