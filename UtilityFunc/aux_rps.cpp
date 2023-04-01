#include <cmath>
#include <vector>
#include "aux.h"
#include <float.h> //for DBL_MAX
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
int* AllocInt(int size)
{
  int* thearray = (int*) malloc(sizeof(int)*size);
  if(thearray==NULL)
  {
    printf("Allocation bombed.\n");
    printf(" size = %3.2f bits (%3.2f MB)\n",(double)size*sizeof(int), sizeof(int)*(double)size / (1024.0*1024.0));
    exit(1);
  }
  int i;
  for (i = 0; i < size; i++){
    thearray[i] = 0;
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
  
  //rps: jobtypes
  //1: A*B
  //2: At*B
  //3: A*Bt
  //4: At*Bt

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
//RPS NOTE:
//This routine returns eigenvectors as a matrix.
//The ROWS are the eigenvectors!
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
//Error-vector-based version from Shepard
//Mol Phys 105 2839 (2007) - eqns 21-23
//Passes error vector, rather than B matrix 
//(thereby eliminating square of condition number)
//Our error vector is a little more complicated
//Rather than [E1 E2 E3 ... En],
//we have [ [E1(q1) E1(q2) ... E1(qN)] [E2(q1) E2(q2) ... E2(qN)] ... [En(q1) En(q2) ... En(qN)] ]
//DIIS coefficients multiply entire sub-[]s.
void linsolverE(double* E, double* B, int N)
{
  //Error vector incoming has dimensions
  // (DIIS subspace size) x (nmodes) x (nbasis*nbasis)
  //Coefficient vector has dimensions
  // (DIIS subspace size)

  //Step 1: Extract last error entry (e_n) 
  //        and build Etilde and ctilde     [eqn 22]
    double* RHS = AllocDouble(N*N);

    /*
    //debug test (9x3 LLS)
    printf("rps debug: dgels test -------------------------------\n");
    double* A = AllocDouble(9*3);
    double* Acopy = AllocDouble(9*3);
    A[0] = 1.000000000000000;   A[1] = 1.361111111111111;   A[2] = 1.911111111111111;
    A[3] = 0.500000000000000;   A[4] = 0.750000000000000;   A[5] = 1.061805555555556;
    A[6] = 0.333333333333333;   A[7] = 0.525000000000000;   A[8] = 0.746203703703704;
    A[9] = 0.500000000000000;   A[10] = 0.750000000000000;   A[11] = 1.061805555555556;
    A[12] = 0.333333333333333;   A[13] = 0.423611111111111;   A[14] = 0.591203703703704;
    A[15] = 0.250000000000000;   A[16] = 0.300000000000000;   A[17] = 0.415902777777778;
    A[18] = 0.333333333333333;   A[19] = 0.525000000000000;   A[20] = 0.746203703703704;
    A[21] = 0.250000000000000;   A[22] = 0.300000000000000;   A[23] = 0.415902777777778;
    A[24] = 0.200000000000000;   A[25] = 0.213611111111111;   A[26] = 0.292722222222222;

    double* b = AllocDouble(9);
    double* bcopy = AllocDouble(9);
    b[0] = 2.690748456790124;
    b[1] = 1.496041666666667;
    b[2] = 1.051729166666667;
    b[3] = 1.496041666666667;
    b[4] = 0.831946373456790;
    b[5] = 0.584916666666667;
    b[6] = 1.051729166666667;
    b[7] = 0.584916666666667;
    b[8] = 0.411254706790123;

    for(int i=0 ; i<9 ; i++)
      bcopy[i] = b[i];
    for(int i=0 ; i<9*3 ; i++)
      Acopy[i] = A[i];

    char trans = 'N';
    int  M     = 9;
         N     = 3;
    int  nrhs  = 1;
    int  lda   = M;
    int  ldb   = M;
    int  lwork = 4*M*N;
    double* work = AllocDouble(lwork);
    int info;
    dgels(&trans,&M,&N,&nrhs,A,&lda,b,&ldb,work,&lwork,&info);

    printf("rps debug: A:\n");
    for(int i=0 ; i<9 ; i++) {
      for(int j=0 ; j<3 ; j++) {
        printf(" % -16.15e",Acopy[i*3+j]);
      }
      printf("\n");
    }
    printf("rps debug: b:\n");
    for(int i=0 ; i<9 ; i++)
      printf(" % -16.15e\n",bcopy[i]);
    printf("rps debug: info = %i\n",info);
    printf("rps debug: bsoln:\n");
    for(int i=0 ; i<9 ; i++)
      printf("% -16.15e\n",b[i]);

    printf("rps debug: A=QR::\n");
    for(int i=0 ; i<9 ; i++) {
      for(int j=0 ; j<3 ; j++) {
        printf(" % -16.15e",A[i*3+j]);
      }
      printf("\n");
    }

    //Check solution (1st 3 elts of b are solns)
      double* test = AllocDouble(9*1);
      for(int i=0 ; i<9 ; i++) {
        for(int j=0 ; j<1 ; j++) {
          for(int k=0 ; k<3 ; k++) {
            test[i*1+j] += Acopy[i*3+k]*b[k*1+j];
          }
        }
      }
      for(int i=0 ; i<9 ; i++)
        printf("rps debug: test[%i] = % -16.15e\n",i,test[i]);
    //end test
    */

    //debug test (3x2 LLS)
    printf("rps debug: dgels test -------------------------------\n");
    double* A = AllocDouble(3*2);
    double* Acopy = AllocDouble(3*2);
    A[0] = 1.0; A[1] = 2.0;
    A[2] = 3.0; A[3] = 4.0;
    A[4] = 5.0; A[5] = 6.0;

//    A[0] = 1.0; A[1] = 3.0;
//    A[2] = 5.0; A[3] = 2.0;
//    A[4] = 4.0; A[5] = 6.0;

    double* b     = AllocDouble(3);
    double* bcopy = AllocDouble(3);
    b[0] = 1.0;
    b[1] = 0.5;
    b[2] = 0.33;

    for(int i=0 ; i<3 ; i++)
      bcopy[i] = b[i];
    for(int i=0 ; i<3*2 ; i++)
      Acopy[i] = A[i];

    char trans = 'N';
    int  M     = 3;
      N     = 2;
    int  nrhs  = 1;
    int  lda   = M;
    int  ldb   = M;
    int  lwork = 4*M*N;
    //double wkopt;
    double* work;
    int info;
    work = AllocDouble(lwork);
    dgels(&trans,&M,&N,&nrhs,A,&lda,b,&ldb,work,&lwork,&info);
/*    int* jpvt= AllocInt(N);
    for(int i=0 ; i<N; i++)
      jpvt[i] = 0;
    double rcond;
    int    rank; 
    dgelsy(&M,&N,&nrhs,A,&lda,b,&ldb,jpvt,&rcond,&rank,work,&lwork,&info);
*/
    printf("rps debug: A:\n");
    for(int i=0 ; i<3 ; i++) {
      for(int j=0 ; j<2 ; j++) {
        printf(" % -16.15e",Acopy[j*3+i]);
      }
      printf("\n");
    }
    printf("rps debug: b:\n");
    for(int i=0 ; i<3 ; i++)
      printf(" % -16.15e\n",bcopy[i]);
    printf("rps debug: info = %i\n",info);
    printf("rps debug: bsoln:\n");
    for(int i=0 ; i<3 ; i++)
      printf("% -16.15e\n",b[i]);

    printf("rps debug: A=QR::\n");
    for(int i=0 ; i<3 ; i++) {
      for(int j=0 ; j<2 ; j++) {
        printf(" % -16.15e",A[i*2+j]);
      }
      printf("\n");
    }

    //Check solution (1st 2 elts of b are solns)
      printf("A*x:\n");
      double* test = AllocDouble(3*1);
      for(int i=0 ; i<3 ; i++) {
        for(int j=0 ; j<1 ; j++) {
          for(int k=0 ; k<2 ; k++) {
            test[i*1+j] += Acopy[k*3+i]*b[k*1+j];
          }
        }
      }
      for(int i=0 ; i<3 ; i++)
        printf(" % -16.15e\n",test[i]);
      printf("A*x-b:\n");
      for(int i=0 ; i<3 ; i++) 
        test[i] -= bcopy[i];
      for(int i=0 ; i<3 ; i++)
        printf(" % -16.15e\n",test[i]);
    //end test
 
    /*
    //debug test (3x3 LLS)
    printf("rps debug: dgels test -------------------------------\n");
    double* A = AllocDouble(3*3);
    double* Acopy = AllocDouble(3*3);
    A[0] = 1.000000000000000;   A[1] = 0.500000000000000;   A[2] = 0.333333333333333;
    A[3] = 0.500000000000000;   A[4] = 0.333333333333333;   A[5] = 0.250000000000000;
    A[6] = 0.333333333333333;   A[7] = 0.250000000000000;   A[8] = 0.200000000000000;

    double* b = AllocDouble(9);
    double* bcopy = AllocDouble(9);
    b[0] = 2.0;
    b[1] = 4.0;
    b[2] = 6.0;

    for(int i=0 ; i<3 ; i++)
      bcopy[i] = b[i];
    for(int i=0 ; i<3*3 ; i++)
      Acopy[i] = A[i];

    char trans = 'N';
    int  M     = 3;
         N     = 3;
    int  nrhs  = 1;
    int  lda   = M;
    int  ldb   = M;
    int  lwork = 2*M*N;
    double* work = AllocDouble(lwork);
    int info;
    dgels(&trans,&M,&N,&nrhs,A,&lda,b,&ldb,work,&lwork,&info);

    printf("rps debug: info = %i\n",info);
    for(int i=0 ; i<3 ; i++)
      printf("rps debug: b[%i] = % -16.15e\n",i,b[i]);

    //Check solution (1st 3 elts of b are solns)
      double* test = AllocDouble(3*1);
      for(int i=0 ; i<3 ; i++) {
        for(int j=0 ; j<1 ; j++) {
          for(int k=0 ; k<3 ; k++) {
            test[i*1+j] += Acopy[i*3+k]*b[k*1+j];
          }
        }
      }
      for(int i=0 ; i<3 ; i++)
        printf("rps debug: test[%i] = % -16.15e\n",i,test[i]);
    //end test
    */

exit(1);

  //Step 2: Solve linear least squares      [eqn 23]
  //dgels() ...or dgelss() for svd
}
//=============================================
void linsolver(double* A, double* B, int N)
{
        //original version:
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

        //rps added:
          //1. Pulay damping
          //  double fac = 0.02;
          //  printf("rps debug: Damping DIIS equations fac = % -3.2f .\n",fac);
          //  for(int i=0 ; i<N-1 ; i++) {
          //    A[i*N+i] *= 1.0+fac;
          //  }

          //2. Normalize to Bii to avoid underflow
          //  printf("rps debug: Normalizing DIIS equations.\n");
          //  double Bii;
          //  for(int i=0 ; i<N-1 ; i++) {
          //    Bii = A[i*N+i];
          //    if(fabs(Bii) > 1.0e-30) {     //avoid div by zero
          //      for(int j=0 ; j<N ; j++) {
          //        //A[i*N+j] /= Bii;
          //        A[i*N+j] /= 1.0;
          //      }
          //    }
          //  }
        //end rps

        dgesvx_(&fact,&trans,&N,&NRS,A,&LDA,AF,&LDA,IPIV,&equed,&R,&C,B,&LDB,X,&LDB,&RCond,Ferr,Berr,Work,IWork,&info);
        printf("rps debug: Rcond = % -10.16e  Info = %i\n",RCond,info); 
        //rps test
        //  printf("rps debug: damping big coefficients\n");
        //  double big = 10.0;
        //  double fac = 0.75;
        //  for(int i=0 ; i<N ; i++) {
        //    if(fabs(X[i]) > big) {
        //      X[i] *= fac;
        //      printf("Scaling X[%i] to % -3.2f\n",i,X[i]);
        //    }
        //  }
        //  double norm = 0.0;
        //  for(int i=0 ; i<N ; i++)
        //    norm += X[i];
        //  for(int i=0 ; i<N ; i++)
        //    X[i] /= norm;
        //end rps
        for(i=0; i<N; i++) B[i] = X[i];

        //////////////////////////////////////////////////////////
        //Shepard E-based version
        //Will need actual error vector instead of B matrix here
        
        //////////////////////////////////////////////////////////

        /*
        //Sellers C2-DIIS version

          //First, isolate Pulay B matrix:
            double* subB = AllocDouble((N-1)*(N-1));
            double* subBcopy = AllocDouble((N-1)*(N-1));
            for(int i=0 ; i<N-1 ; i++) 
              for(int j=0 ; j<N-1 ; j++) 
                subB[i*(N-1)+j] = A[i*N+j];
            printf("Original DIIS:\n");
            for(int i=0 ; i<N ; i++) {
              for(int j=0 ; j<N ; j++)
                printf(" % -10.8e",A[i*N+j]);
              printf("\n");
            }
            printf("subDIIS:\n");
            for(int i=0 ; i<N-1 ; i++) {
              for(int j=0 ; j<N-1 ; j++)
                printf(" % -16.15e",subB[i*(N-1)+j]);
              printf("\n");
            }

          //Then diagonalize
          N--;
          double* val = AllocDouble(N);
          diagonalize(subB,val,N); //evecs as ROWS
          printf("Eigenvalues of B:\n");
          for(int i=0 ; i<N ; i++)
            printf("% -16.15e\n",val[i]);
          //printf("Eigenvectors of B (as rows):\n");
          //for(int i=0 ; i<N ; i++) {
          //  for(int j=0 ; j<N ; j++) 
          //    printf(" % -16.15e",subB[i*N+j]);
          //  printf("\n");
          //}

          //1-norm the eigenvectors
          double norm = 0.0;
          for(int i=0 ; i<N ; i++) {
            norm = 0.0;
            for(int j=0 ; j<N ; j++) {
              norm += subB[i*N+j];
            }
            //printf("%i norm = % -8.4f\n",i,norm);
            for(int j=0 ; j<N ; j++) {
              subB[i*N+j] /= norm;
            }
          }
          printf("Normalized eigenvectors of B (as rows):\n");
          for(int i=0 ; i<N ; i++) {
            for(int j=0 ; j<N ; j++) 
              printf(" % -16.15e",subB[i*N+j]);
            printf("\n");
          }

          //Compute estimated errors from each eigenvector
            double* err = AllocDouble(N);
            double* temp = AllocDouble(N);
            for(int k=0 ; k<N ; k++) { //loop over solutions
              err[k] = 0.0;
              for(int i=0 ; i<N ; i++)
                temp[i] = subB[k*N+i];
              for(int i=0 ; i<N ; i++)
                for(int j=0 ; j<N ; j++)
                  err[k] += temp[i]*subB[i*N+j]*temp[j];
              printf("err %i = % -16.15e\n",k,err[k]);
            }

          //Add: Cap on max elts

          //Save best coefficients
            double minerr = DBL_MAX;
            int    minerr_idx = -1;
            for(int k=0 ; k<N ; k++) {
              if(fabs(err[k]) < fabs(minerr)) {
                minerr = err[k];
                minerr_idx = k;
              }
            }
            printf("min err          = % -16.15e\n",minerr);
            printf("min err solution = %i\n",minerr_idx);
            for(int i=0 ; i<N ; i++)
              B[i] = subB[minerr_idx*N+i];
        */

        //////////////////////////////////////////////////////////

        /*
        //SVD version: [KEEP]
          //First, isolate Pulay B matrix:
            double* subB = AllocDouble((N-1)*(N-1));
            double* subBcopy = AllocDouble((N-1)*(N-1));
            for(int i=0 ; i<N-1 ; i++) 
              for(int j=0 ; j<N-1 ; j++) 
                subB[i*(N-1)+j] = A[i*N+j];
          printf("Original DIIS:\n");
          for(int i=0 ; i<N ; i++) {
            for(int j=0 ; j<N ; j++)
              printf(" % -10.8e",A[i*N+j]);
            printf("\n");
          }
          printf("subDIIS:\n");
          for(int i=0 ; i<N-1 ; i++) {
            for(int j=0 ; j<N-1 ; j++)
              printf(" % -16.15e",subB[i*(N-1)+j]);
            printf("\n");
          }

          //Then SVD
          N--;

            //rps added: normalize to Bii to avoid underflow
            //  printf("rps debug: Normalizing DIIS equations.\n");
            //  double Bii;
            //  for(int i=0 ; i<N ; i++) {
            //    Bii = subB[i*N+i];
            //    if(fabs(Bii) > 1.0e-30) {     //avoid div by zero
            //      for(int j=0 ; j<N ; j++) {
            //        subB[i*N+j] /= Bii;
            //      }
            //    }
            //  }

            //  double fac = 0.02;
            //  printf("rps debug: Damping DIIS equations fac = % -3.2f.\n",fac);
            //  for(int i=0 ; i<N ; i++) 
            //    subB[i*N+i] *= 1.0 + fac;
            //end rps
            
          for(int i=0 ; i<N*N ; i++)
            subBcopy[i] = subB[i];
          char   jobu      = 'A'; //all vectors
          char   jobvt     = 'A'; //all vectors
          int    M         = N;
          int    LDA       = N;
          double S[N];
          double U[N*N];
          double Vt[N*N];
          int    LDU       = N;
          int    LDVt      = N;
          double Work[5*N];
          int    LWork     = 5*N;
          int    info;

          dgesvd_(&jobu,&jobvt,&M,&N,subB,&LDA,S,U,&LDU,Vt,&LDVt,Work,&LWork,&info);
          printf("rps debug: SVD condition number: % -8.4e\n",S[0]/S[N-1]);
            printf("U:\n");
            for(int i=0 ; i<N ; i++) {
              for(int j=0 ; j<N ; j++)
                printf(" % -10.8e",U[i*N+j]);
              printf("\n");
            }
          for(int i=0 ; i<N ; i++)
            printf("rps debug: singval[%i] = % -16.15e\n",i,S[i]);
            printf("Vt:\n");
            for(int i=0 ; i<N ; i++) {
              for(int j=0 ; j<N ; j++)
                printf(" % -10.8e",Vt[i*N+j]);
              printf("\n");
            }

          //Build pseudoinverse
            double cutoff = 1.0e-11; //cutoff is effectively 1/condition number.  Multiplies sigma_max
            double* pInv = AllocDouble(N*N);
            double* Winv = AllocDouble(N*N);
            for(int i=0 ; i<N ; i++) {
              for(int j=0 ; j<N ; j++) {
                if (i==j) {
                  if(fabs(S[i]) > S[0]*cutoff)
                    Winv[i*N+j] = 1/S[i];
                  else
                  {
                    printf("Skipping singval = % -16.15e\n",S[i]);
                    Winv[i*N+j] = 0.0;
                  }
                }
                else
                  Winv[i*N+j] = 0.0;
              }
            }
            double* temp = AllocDouble(N*N); //Winv * Ut
            ABmult(temp,Winv,U,N,N,N,N,N,N,3);
            printf("temp:\n");
            for(int i=0 ; i<N ; i++) {
              for(int j=0 ; j<N ; j++)
                printf(" % -10.8e",temp[i*N+j]);
              printf("\n");
            }
            ABmult(pInv,Vt,temp,N,N,N,N,N,N,2);
            printf("pInv:\n");
            for(int i=0 ; i<N ; i++) {
              for(int j=0 ; j<N ; j++)
                printf(" % -10.8e",pInv[i*N+j]);
              printf("\n");
            }
            //Test inverse:
            ABmult(temp,pInv,subBcopy,N,N,N,N,N,N,1);
            printf("A^-1*A:\n");
            for(int i=0 ; i<N ; i++) {
              for(int j=0 ; j<N ; j++)
                printf(" % -10.8e",temp[i*N+j]);
              printf("\n");
            }
          //Solve Ax=b --> x=(A^-1)b via pseudoinverse
          //Use b=1 to match Q-Chem [not sure why yet...]
            double* b = AllocDouble(N);
            for(int i=0 ; i<N ; i++)
              b[i] = 1.0;
            double* soln = AllocDouble(N);
            for(int i=0 ; i<N ; i++) 
              soln[i] = 0.0;
            for(int i=0 ; i<N ; i++) 
              for(int j=0 ; j<N ; j++) 
                soln[i] += pInv[i*N+j]*b[j];
            printf("A^-1*b:\n");
            for(int i=0 ; i<N ; i++) 
              printf("% -16.15e\n",soln[i]);
            printf("1-normalized soln:\n");
            double norm1 = 0.0;
            for(int i=0 ; i<N ; i++) 
              norm1 += soln[i];
              //norm += pow(soln[i],2);
            for(int i=0 ; i<N ; i++) 
              printf("% -16.15e\n",soln[i]/norm1);

            for(int i=0 ; i<N ; i++) 
              B[i] = soln[i]/norm1;
        */

        /*
        //temp test first:
        //Test works...shows reduction in RCond with reduction in subspace
        double* Acopy = AllocDouble(N*N);
        for(int i=0 ; i<N*N ; i++) 
          Acopy[i] = A[i];
        double* Xcopy = AllocDouble(N);
        for(int s=0 ; s<N-1 ; s++) {

          for(int i=0 ; i<N*N ; i++)
            A[i] = 0.0;

          //Build sub-DIIS matrix from original
          //B part first
          for(int i=0 ; i<(N-s)-1 ; i++) 
            for(int j=0 ; j<(N-s)-1 ; j++) 
              A[i*(N-s)+j] = Acopy[i*N+j];
          //Then the 1's
          for(int i=0 ; i<(N-s)-1 ; i++) {
            //Last col
            A[i*(N-s)+(N-s-1)] = -1.0;
            //Last row
            A[(N-s-1)*(N-s)+i] = -1.0;
          }

          printf("Original DIIS matrix:\n");
          for(int i=0 ; i<N ; i++) {
            for(int j=0 ; j<N ; j++) 
              printf(" % -10.8e",Acopy[i*N+j]);
            printf("\n");
          }
          printf("Reduced DIIS matrix:\n");
          for(int i=0 ; i<N-s ; i++) {
            for(int j=0 ; j<N-s ; j++) 
              printf(" % -10.8e",A[i*(N-s)+j]);
            printf("\n");
          }

          int i;
          int NRS = 1;
          int LDA = N-s;
          int LDB = N-s;
          double Berr[1];
          double X[N-s];
          int IPIV[N-s];
          int IWork[N-s];
          double Ferr[1];
          double Work[4*(N-s)];
          double AF[LDA*(N-s)];
          double RCond, R, C;
          char fact = 'N';
          char trans= 'T';
          char equed = 'X';
          int info;
          int Ntemp = N-s;
          dgesvx_(&fact,&trans,&Ntemp,&NRS,A,&LDA,AF,&LDA,IPIV,&equed,&R,&C,B,&LDB,X,&LDB,&RCond,Ferr,Berr,Work,IWork,&info);
          printf("rps debug: DIIS test subspace loses %i RCond = % -4.2e\n",s,RCond);

          if(s==0) {
            for(int i=0 ; i<N; i++)
              Xcopy[i] = X[i];
          }
        }
        for(int i=0; i<N; i++) B[i] = Xcopy[i];
        */


        /*
        //New: Solve DIIS equations
        //If condition number is problematic, reduce subspace size and try again
        //Use empirical choice of inverse condition number:
        double cutoff = 1e-13;
        int ok = 0;
        int s  = 0;
        double* Acopy = AllocDouble(N*N);
        for(int i=0 ; i<N*N ; i++) 
          Acopy[i] = A[i];

        while (ok==0) {

          //Build (possibly adjusted) DIIS matrix
          double* Anew = AllocDouble((N-s)*(N-s));
          for(int i=0 ; i<(N-s)*(N-s) ; i++)
            Anew[i] = 0.0;
          //Build sub-DIIS matrix from original
          //B part first
          //Lower right:
          for(int i=0 ; i<(N-s)-1 ; i++) 
            for(int j=0 ; j<(N-s)-1 ; j++) 
              Anew[i*(N-s)+j] = Acopy[(i+s)*N+(j+s)];
          //Then the 1's
          for(int i=0 ; i<(N-s)-1 ; i++) {
            //Last col
            Anew[i*(N-s)+(N-s-1)] = -1.0;
            //Last row
            Anew[(N-s-1)*(N-s)+i] = -1.0;
          }

            printf("Original DIIS matrix:\n");
            for(int i=0 ; i<N ; i++) {
              for(int j=0 ; j<N ; j++) 
                printf(" % -10.8e",Acopy[i*N+j]);
              printf("\n");
            }
            printf("Reduced DIIS matrix:\n");
            for(int i=0 ; i<N-s ; i++) {
              for(int j=0 ; j<N-s ; j++) 
                printf(" % -10.8e",Anew[i*(N-s)+j]);
              printf("\n");
            }


          //Build (possibly adjusted) RHS
          double* Bnew = AllocDouble(N-s);
          for(int i=0 ; i<(N-s-1) ; i++)
            Bnew[i] = 0.0;
          Bnew[N-s-1] = -1.0;
          printf("rps debug: Bnew:\n");
          for(int i=0 ; i<(N-s) ; i++)
            printf("% -16.15e\n",Bnew[i]);


          //Solve DIIS equations in this subspace
          int i;
          int NRS = 1;
          int LDA = N-s;
          int LDB = N-s;
          double Berr[1];
          double X[N-s];
          int IPIV[N-s];
          int IWork[N-s];
          double Ferr[1];
          double Work[4*(N-s)];
          double AF[LDA*(N-s)];
          double RCond, R, C;
          char fact = 'N';
          char trans= 'T';
          char equed = 'X';
          int info;
          int Ntemp = N-s;
          dgesvx_(&fact,&trans,&Ntemp,&NRS,Anew,&LDA,AF,&LDA,IPIV,&equed,&R,&C,Bnew,&LDB,X,&LDB,&RCond,Ferr,Berr,Work,IWork,&info);
          printf("rps debug: RCond = % -4.2e\n",RCond);
          for(int i=0 ; i<(N-s) ; i++)
            printf("rps debug: X[%i] = % -4.2e\n",i,X[i]);

          if(RCond > cutoff) {
            ok=1;
            for(int i=0; i<N; i++) {
              if(i<s)
                B[i] = 0.0;
              else
                B[i] = X[i-s];
            }
          }
          else {
            s++;
            printf("RCond too small ( % -4.2e)\n",RCond);
            printf("Adjusting subspace size to %i\n",N-s-1);
          }

          free(Anew);
          free(Bnew);
        }
        */
}

void printmat(double* mat, int o, int n, int m, double scale)
{
int i,j,k;
for(k=0 ; k<o ; k++) {
for(i=0; i<n; i++){
  for(j=0; j<m; j++){
    printf("% -10.6e  ",mat[k*n*m+i*m+j]*scale);
  }
 printf("\n");
 }
printf("\n");
}
printf("\n");
}

int tupleIndex(int& counter, int dim, int recursionLevel, std::vector<int>& indices, std::vector<int>& targetIndices, int nModes) {
  if(recursionLevel == dim-1)
    while(indices[recursionLevel] < nModes-1) {
      //Check if indices match
      if(indices == targetIndices)
        return counter;
      //Increment if not
      indices[recursionLevel]++;
      counter++;  
  } else {
    while(indices[recursionLevel] < nModes-indices.size()+recursionLevel) {
      tupleIndex(counter,dim,recursionLevel+1,indices,targetIndices,nModes);
      if(indices == targetIndices)
        return counter;
      indices[recursionLevel]++;
      for(int i=recursionLevel+1; i<dim ; i++) 
        indices[i] = indices[i-1]+1;
      counter++;
    }
  } 
  return counter;
}

int tupleIndexDriver(std::vector<int>& targetIndices, int nModes) {
  for(int i=0 ; i<targetIndices.size() ; i++) {
    if(targetIndices[i] > nModes-1)
      return -1;
  }
  int dim = targetIndices.size();
  std::vector<int> indices(dim); 
  for(int i=0 ; i<dim ; i++)
    indices[i] = i;
  int counter = 0;
  
  return tupleIndex(counter,dim,0,indices,targetIndices,nModes);
}
