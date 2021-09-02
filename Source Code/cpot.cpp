#include "cpot.h"
#include<iostream>
#include<fstream>
#include<vector>
#include<string.h>
#include<cmath>
#include <sstream>

using namespace std;

//cpot::cpot(double* pots, int degrees, int grid_sizes, int dofs){
//double* pot = pots;
cpot::cpot(int grid_sizes) { 
grid_size = grid_sizes;
//degree = degrees;
//dof = dofs;
}

cpot::~cpot() {}
void cpot::get_coupling(vector<double>& pot, int degree) {
if( degree ==6)
get6d(pot);
else if( degree ==5)
get5d(pot);
else if( degree == 4)
get4d(pot);
else if (degree == 3)
//get3d(pot,grid_size,dof);
get3d(pot);
else
////get2d(pot,grid_size,dof);
get2d(pot);
}
//


//double cpot::get2d(double* pot2,int grid_size,int dof){
void cpot::get2d(vector<double>& pot2) {
  const char*hello = "billy" ;
  //int len = (fact(dof)/(fact(dof-2)*fact(2)))*(int)(pow(grid_size,2));
  int len = (int)(pow(grid_size,2));
  //printf("%i",len);

  int sub = (grid_size-1)/2;
  double equilibrium = pot2[(((int)pow(grid_size,2)-1)/2)];
//  printf("\n%f\n",equilibrium);
//  for(int i = 0; i < len; i++){
//    pot2[i] = pot2[i] - equilibrium ;
//  }

  int n = 2;
  vector<int> index(len);
  int** stored_values = new int*[len];
  for(int i=0 ; i<len ; i++) {
    stored_values[i] = new int[2];
    //printf("\n%f\n",pot2[i]);
  }
  int count = -1;
  int factor = 0;
  vector<double> twoD(len);
  for(int i = 0;i<len;i++){twoD[i] = pot2[i];}

//for(int var1 = 1; var1 < 2; var1++){
//	for(int var2 = var1+1; var2 <3; var2++){
		for(int var3 =-sub; var3 <sub+1;var3++){
			for(int var4 =-sub; var4 <sub+1;var4++){
				count = count +1;
//				cout << count << '\n';
				int term1 = var4;
				int term2 = var3 * grid_size;
//				int term3 = (var2-2)*(int)pow(grid_size,2);
//				int term4 = factor*(int)pow(grid_size,2);
				stored_values[count][0] = term1;
				stored_values[count][1] = term2;
		//		stored_values[count][2] = term3;
		//		stored_values[count][3] = term4;
					
				for(int in = 0; in <2; in++){
				//printf("\n%i\n",stored_values[count][in]);
				  index[count] +=  stored_values[count][in];
				//printf("\n%i\n",index[count]);
		    }//for in
			}//for var4
		}//for var3
//	}
//    factor = factor + (dof - (var1+1));
//}
//for(int i=0;i<121;i++){
//printf("\n%f\n",twoD[i]);
//}


  int** dummy = new int*[len];
  for(int i=0 ; i<len ; i++) {
    dummy[i] = new int[2];
  }
  for(int q = 0;q < len; q++){
	  for(int w=0; w<2;w++){
	    dummy[q][w] = stored_values[q][w];
    //	printf("\n%i\n",stored_values[q][w]);
    }
  }

  int corr = -index[0];
//  printf("\n%i\n",corr);
//  for(int u = 0; u <len;u++){
  //    printf("\n%f\n",pot2[u]);
//  }


  for(int d = 0; d<2; d++){
	  for(int l = 0;l<len;l++){
		  int new_index = 0;
		  for(int o = 0;o<2;o++){
		    stored_values[l][o] = dummy[l][o];
		  }
		  stored_values[l][d] = 0;
		  for( int k = 0; k<2; k++){
			  new_index = new_index +  stored_values[l][k];
				//printf("\n%i\n",new_index);
		  }
	    new_index = new_index + corr;
//	printf("\n%i\n",new_index);
	    pot2[(index[l])+corr]  = pot2[(index[l])+corr]-twoD[(new_index)];
    }//for l
  }//for d

  for(int i=0 ; i<len ; i++) {delete[] stored_values[i];}
  delete[] stored_values;
  for(int i=0 ; i<len ; i++) {delete[] dummy[i];}
  delete[] dummy;

//int p ;
//for(p=0;p<len;p++) printf("\n%f\n",twoD[p]);
//return *twoD;*/
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//double cpot::get3d(double* pot3,int grid_size,int dof){
void cpot::get3d(vector<double>& pot3) {
  const char*hello = "billy" ;
  for(int k =0;k<1331;k++){
//printf("\n%f\n",pot3[k]);
  }
//printf("\n%i\n",kip);
  int len = (int)(pow(grid_size,3));
  int sub = (grid_size-1)/2;
//  printf("%i",len);
  double equilibrium = pot3[(((int)pow(grid_size,3)-1)/2)];
//  printf("\n%f\n",equilibrium);

  for(int i = 0; i < len; i++){
    pot3[i] = pot3[i] - equilibrium ;
  }
  int n = 3;
  vector<int> index(len);
  int** stored_values = new int*[len];
  for(int i=0 ; i<len ; i++) {
    stored_values[i] = new int[3];
  }

  int count = -1;
  int factor1 = 0;
  int factor2 = 0;
  vector<double> threeD(len);
	for(int var4 = -sub; var4 < sub+1;var4++){
	  for(int var5 = -sub; var5 < sub+1;var5++){
		  for(int var6 = -sub; var6<sub+1;var6++){
			  count = count +1;
//			cout << count << '\n';
				int term1 = var6;
				int term2 = var5 * (grid_size);
				int term3 = var4 * ((int)pow(grid_size,2));
				//int term4 = (var3-3)*((grid_size)^3);
				//int term5 = (factor2) * ((grid_size)^3);
				//int term6 = (factor1) * ((grid_size)^3);
				stored_values[count][0] = term1;
				stored_values[count][1] = term2;
				stored_values[count][2] = term3;

				for(int in = 0; in <3; in++){
          index[count] = index[count] + stored_values[count][in];
				}
			}
		}
  }

  int** dummy = new int*[len];
  for(int i=0 ; i<len ; i++) {
    dummy[i] = new int[3];
//  printf("\n%i\n",index[i]);
}
  for(int q = 0;q < len; q++){
    for(int w=0; w<3;w++){
      dummy[q][w] = stored_values[q][w];
	//printf("\n%i\n",dummy[q][w]);
    }
  }

  for(int i = 0;i<len;i++){
	  threeD[i] = pot3[i];
//	printf("\n%f\n",pot3[i]);
  }

  int corr = -index[0];
/*for(int d = 0; d<3; d++){
        for(int l = 0;l<len;l++){
                int new_index = 0;
                for(int o = 0;o<3;o++){
                stored_values[l][o] = dummy[l][o];
                }
                stored_values[l][d] = 0;
                for( int k = 0; k<3; k++){
  //                   new_index +=  stored_values[l][k];
//		}
	new_index = new_index + corr;
//	printf("\n%i\n",new_index);
//        threeD[(index[l])+corr]  = threeD[(index[l])+corr]-pot3[(new_index)];
        }
	if(d==2){*/
//	for(int i = 0; i <len; i++){
	//	pot3[i] = threeD[i];
	//}
//	for(int d = 0; d<2;d++){
//		for(int dd = d+1; dd <3;dd++){
//
                for(int d = 0; d < 3; d++){
			for(int m=0;m<len;m++){
				int new_index =0;
				for(int o = 0; o<3;o++){
				  stored_values[m][o] = 0;
					}
				
				  stored_values[m][d] = dummy[m][d];
				for(int k = 0; k<3;k++)
				  new_index = new_index+ stored_values[m][k];
				new_index = new_index + corr;
//				 printf("\n%i\n",new_index);	
			  pot3[index[m]+corr]  = pot3[(index[m])+corr]-threeD[(new_index)];	
      }//for m
		  if(d==2){
			  for(int i = 0; i<len;i++)
			    threeD[i] = pot3[i];
			  //printf(hello);
			  for(int d=0; d<3; d++){
          for(int l = 0;l<len;l++){
            int new_index = 0;
            for(int o = 0;o<3;o++)
              stored_values[l][o] = dummy[l][o];
            stored_values[l][d] = 0;
            for(int k = 0; k<3; k++)
              new_index +=  stored_values[l][k];
            new_index = new_index + corr;
                        //    printf("\n%i\n",new_index);
  		      pot3[(index[l])+corr]  = pot3[(index[l])+corr]-threeD[(new_index)];               
          }//for l
        }//for d
      }//if d==1 d==2
    }//for dd
  //}//for d */
//}
  for(int i=0 ; i<len ; i++) {delete[] stored_values[i];}
  delete[] stored_values;
  for(int i=0 ; i<len ; i++) {delete[] dummy[i];}
  delete[] dummy;
//for(int q = 0; q <len;q++){printf("\n%f %i\n",pot3[q],q);}
//Run Emily Code to convert to 2d pots
//Call get2D
//retrun threeD
}//end method

////////////////////////////////////////////////////////////////////////////////////////////////

void cpot::get4d(vector<double>& pot4){
const char*hello = "billy" ;

int len = (int)(pow(grid_size,4));
//printf("\n%d\n",len);
int sub = (grid_size-1)/2 ;

double equilibrium = pot4[(((int)pow(grid_size,4)-1)/2)];

  for(int i = 0; i < len; i++){
    pot4[i] = pot4[i] - equilibrium ;
  }
//for(int q = 0; q <len;q++){printf("\n%f %i\n",pot4[q],q);}
vector<int> index(len);
int** stored_values = new int*[len];
for(int i=0; i<len;i++){
stored_values[i] = new int[4];
}
int count = -1;
vector<double> fourD(len);
for(int var5 = -sub; var5 <sub+1; var5++){
	for(int var6 = -sub; var6 <sub+1; var6++){
		for(int var7 = -sub; var7 <sub+1;var7++){
			for(int var8 = -sub; var8 < sub+1; var8++){
			count = count + 1;
			int term1 = var8;
			int term2 = var7 * (grid_size);
			int term3 = var6 * ((int)pow(grid_size,2));
			int term4 = var5 * ((int)pow(grid_size,3));
			stored_values[count][0] = term1;
			stored_values[count][1] = term2;
			stored_values[count][2] = term3;
			stored_values[count][3] = term4;

			for(int in = 0; in<4; in++){
			index[count] = index[count] + stored_values[count][in];
			}
		}
	}	
}
}

int** dummy = new int*[len];
for(int i=0; i<len;i++){
dummy[i] = new int[4];
}
for(int q = 0; q < len; q++){
	for(int w = 0; w<4; w++){
		dummy[q][w] = stored_values[q][w];
	}
}
for(int i = 0; i<len; i++){
	fourD[i] = pot4[i];
}

int corr = -index[0];

/*

for(int d = 0; d<2; d++){
	for(int dd = d+1;dd<3;dd++){
		for(int ddd = dd+1;ddd<4;ddd++){
			for(int m = 0; m<len;m++){
				int new_index = 0;
				for(int o = 0; o<4;o++)
					stored_values[m][o] = dummy[m][o];

					stored_values[m][d] = 0;
					stored_values[m][dd] = 0;
					stored_values[m][ddd] = 0;
					for(int k = 0; k<4; k++)
					new_index = new_index + stored_values[m][k];
					
					new_index = new_index + corr;
					pot4[index[m]+corr] = pot4[index[m]+corr]-fourD[new_index];
				} //for m
			if(d==1
			
*/
int* looper = new int[4];
for(int i=0; i<4; i++){
	looper[i] = i;
}
for(int j = 0;j<4;j++){
for(int m=0; m< len; m++){
int new_index = 0;
for(int u = 0; u<4;u++){
	stored_values[m][u]=0;
}
stored_values[m][j] = dummy[m][j];

for(int k = 0; k<4; k++)         
new_index = new_index + stored_values[m][k];
                                        
new_index = new_index + corr;

pot4[index[m]+corr] = pot4[index[m]+corr]-fourD[new_index];
}
if(j==3){
for(int i = 0; i<len;i++)
fourD[i] = pot4[i];

for(int d = 3; d > 0; d--){
	for(int dd = d-1; dd >-1; dd--){
		for(int l = 0; l <len; l++){
			int new_index = 0;
			for(int u = 0; u<4;u++){
				stored_values[l][u] = dummy[l][u];
				}
			stored_values[l][d] = 0;
			stored_values[l][dd] = 0;
			for(int k = 0;k<4;k++){
				new_index = new_index + stored_values[l][k];
			}
				new_index = new_index + corr;
			pot4[index[l]+corr] = pot4[index[l]+corr]-fourD[new_index];

			} //for l
		//	printf("\n%i %i\n",d,dd);
		if(d ==1 && dd==0){
//			printf(hello);
			for(int i = 0; i<len;i++)
				fourD[i] = pot4[i];
			for(int c = 3; c > -1; c--){
				for(int n = 0; n <len; n++){		
					int new_index = 0;
				 	for(int u = 0; u<4;u++){
						stored_values[n][u] = dummy[n][u];
					}
					stored_values[n][c] = 0;
					for(int k = 0 ; k<4;k++){
						new_index = new_index + stored_values[n][k];
						}
					new_index = new_index + corr;
				pot4[index[n]+corr] = pot4[index[n]+corr]-fourD[new_index];
			}//for n
}// for c		
}//if d==1,dd==2
} //for dd
}//for d
}//if j==3

}//for j


  for(int i=0 ; i<len ; i++) {delete[] stored_values[i];}
  delete[] stored_values;
  for(int i=0 ; i<len ; i++) {delete[] dummy[i];}
  delete[] dummy;

//for(int q = 0; q <len;q++){printf("\n%f %i\n",pot4[q],q);}

}//end method


///////////////////////////////////////////////////////////////////////////////

void cpot::get5d(vector<double>& pot5){

int len = (int)(pow(grid_size,5));
//printf("\n%d\n",len);
int sub = (grid_size-1)/2 ;
vector<int> index(len);

double equilibrium = pot5[(((int)pow(grid_size,5)-1)/2)];

  for(int i = 0; i < len; i++){
    pot5[i] = pot5[i] - equilibrium ;
  }

int** stored_values = new int*[len];
for(int i=0; i<len;i++){
stored_values[i] = new int[5];
}


//for(int q = 0; q <len;q++){printf("\n%f %i\n",pot5[q],q);}
int count = -1;
vector<double> fiveD(len);

for(int var5 = -sub; var5 <sub+1; var5++){
        for(int var6 = -sub; var6 <sub+1; var6++){
                for(int var7 = -sub; var7 <sub+1;var7++){
                        for(int var8 = -sub; var8 < sub+1; var8++){
                        	for(int var9 = -sub; var9 <sub+1;var9++){
		        	count = count + 1;
                        	int term1 = var9;
                        	int term2 = var8 * (grid_size);
                        	int term3 = var7 * ((int)pow(grid_size,2));
                        	int term4 = var6 * ((int)pow(grid_size,3));
				int term5 = var5 * ((int)pow(grid_size,4));
                       		 stored_values[count][0] = term1;	
                       		 stored_values[count][1] = term2;
                       		 stored_values[count][2] = term3;
                      		  stored_values[count][3] = term4;
				stored_values[count][4] = term5;


                        for(int in = 0; in<5; in++){
                        index[count] = index[count] + stored_values[count][in];
                        }
			}
                }
        }
}
}
int** dummy = new int*[len];
for(int i=0; i<len;i++){
dummy[i] = new int[5];
}

for(int q = 0; q<len;q++){
	for(int a = 0; a<5;a++){
		dummy[q][a] = stored_values[q][a];
		fiveD[q] = pot5[q];
}
}

int corr = -index[0];
int* looper = new int[5];
for(int i=0; i<5; i++){
        looper[i] = i;
}
for(int j = 0;j<5;j++){
for(int m=0; m< len; m++){
int new_index = 0;
for(int u = 0; u<5;u++){
        stored_values[m][u]=0;
}
stored_values[m][j] = dummy[m][j];

for(int k = 0; k<5; k++){
new_index = new_index + stored_values[m][k];
}

new_index = new_index + corr;

pot5[index[m]+corr] = pot5[index[m]+corr]-fiveD[new_index];
}//m
//printf("\n%i\n",j);
if(j==4){
	for(int i=0; i <len; i++){
		fiveD[i] = pot5[i];	
	}
	for(int d = 4; d > 1; d--){
		for(int dd = d-1;dd>0;dd--){
			for(int ddd = dd-1;ddd>-1;ddd--){
				for(int l=0;l<len;l++){
					int new_index = 0;
					for(int u = 0; u<5;u++){
					stored_values[l][u]=dummy[l][u];}//u
				
					stored_values[l][d]=0;
					stored_values[l][dd]=0;
					stored_values[l][ddd]=0;
					for(int k =0;k<5;k++){
						new_index = new_index + stored_values[l][k];
						}//k
						new_index = new_index + corr;
						pot5[index[l]+corr] = pot5[index[l]+corr]-fiveD[new_index];
}//l
		//printf("\n%i %i %i\n",d,dd,ddd);
if(d ==2 && dd==1 && ddd==0){

        for(int i=0; i <len; i++){
                fiveD[i] = pot5[i];
        }
	for(int c = 4; c>0; c--){
		for(int cc = c-1;cc>-1;cc--){
			for(int n = 0; n <len; n++) {
				int new_index = 0;
				for(int u = 0; u<5;u++){
					stored_values[n][u] = dummy[n][u];}//u
					stored_values[n][c] = 0;
					stored_values[n][cc] = 0;
					for(int k = 0; k<5; k++){
					new_index = new_index + stored_values[n][k];
					}
					new_index = new_index + corr;
					pot5[index[n]+corr] = pot5[index[n]+corr]-fiveD[new_index];			

}//n
			//printf("\n%i %i\n",c,cc);

if(c==1 && cc==0){
        for(int i=0; i <len; i++){
                fiveD[i] = pot5[i];
        }
	for(int b = 4;b>-1;b--){
		for(int p = 0; p<len;p++){
			int new_index = 0;
			for(int u=0; u<5;u++){	
				stored_values[p][u] = dummy[p][u];}//u
				stored_values[p][b] = 0;
				for(int k = 0; k<5;k++){
				new_index = new_index + stored_values[p][k];
				}
				new_index = new_index + corr;
				pot5[index[p]+corr] = pot5[index[p]+corr]-fiveD[new_index];
}//p
	//	printf("\n%i\n",b);
}//b
}//if c==1,cc==0

}//cc

}//c


	

}//if d =2,dd=1,ddd=0
}//ddd
}//dd



}//d
}//if j==4*/
}//j

  for(int i=0 ; i<len ; i++) {delete[] stored_values[i];}
  delete[] stored_values;
  for(int i=0 ; i<len ; i++) {delete[] dummy[i];}
  delete[] dummy;

//for(int q = 0; q <len;q++){printf("\n%f %i\n",pot5[q],q);}
}//end of method



//////////////////////////////////////////////////////////////////////////////



void cpot::get6d(vector<double>& pot6){
int len = (int)(pow(grid_size,6));
//printf("\n%d\n",len);
int sub = (grid_size-1)/2 ;
vector<int> index(len);

double equilibrium = pot6[(((int)pow(grid_size,6)-1)/2)];

  for(int i = 0; i < len; i++){
    pot6[i] = pot6[i] - equilibrium ;
  }

int** stored_values = new int*[len];
for(int i=0; i<len;i++){
stored_values[i] = new int[6];
}

int count = -1;
vector<double> sixD(len);


for(int var5 = -sub; var5 <sub+1; var5++){
        for(int var6 = -sub; var6 <sub+1; var6++){
                for(int var7 = -sub; var7 <sub+1;var7++){
                        for(int var8 = -sub; var8 < sub+1; var8++){
                                for(int var9 = -sub; var9 <sub+1;var9++){
                               		for(int var10 = -sub; var10 < sub+1; var10++){
						 count = count + 1;
						 int term1 = var10;
						 int term2 = var9 * (grid_size);
						 int term3 = var8 * ((int)pow(grid_size,2));
						 int term4 = var7 * ((int)pow(grid_size,3));
					         int term5 = var6 * ((int)pow(grid_size,4));
						 int term6 = var5 * ((int)pow(grid_size,5));
						  stored_values[count][0] = term1;
                                 stored_values[count][1] = term2;
                                 stored_values[count][2] = term3;
                                  stored_values[count][3] = term4;
                                stored_values[count][4] = term5;
				stored_values[count][5] = term6;
				                        for(int in = 0; in<6; in++){
                        index[count] = index[count] + stored_values[count][in]; 
			}
			}
		}
	}
	}	
	}
	}





int** dummy = new int*[len];
for(int i=0; i<len;i++){
dummy[i] = new int[6];
}



for(int q = 0; q<len;q++){
        for(int a = 0; a<6;a++){
                dummy[q][a] = stored_values[q][a];
                sixD[q] = pot6[q];
}
}

int corr = -index[0];

//printf("\n%d\n",len);

for(int j = 0;j<6;j++){
for(int m=0; m< len; m++){
int new_index = 0;
for(int u = 0; u<6;u++){
        stored_values[m][u]=0;
}
stored_values[m][j] = dummy[m][j];

for(int k = 0; k<6; k++){
new_index = new_index + stored_values[m][k];
}

new_index = new_index + corr;

pot6[index[m]+corr] = pot6[index[m]+corr]-sixD[new_index];
}//m
//printf("\n%i\n",j);
if(j==5){
        for(int i=0; i <len; i++){
                sixD[i] = pot6[i];
        }
       for(int d = 5; d > 2; d--){
                for(int dd = d-1;dd>1;dd--){
                        for(int ddd = dd-1;ddd>-0;ddd--){
				for(int dddd = ddd-1;dddd>-1;dddd--){
					for(int l = 0; l<len;l++){
						int new_index = 0;
						for(int u=0;u<6;u++){
							stored_values[l][u] = dummy[l][u];
							}
						stored_values[l][d] = 0;
						stored_values[l][dd]=0;
						stored_values[l][ddd]=0;
						stored_values[l][dddd]=0;
						for(int k = 0; k<6;k++){
							new_index = new_index + stored_values[l][k];
							}
							new_index = new_index + corr;
							pot6[index[l]+corr] = pot6[index[l]+corr]-sixD[new_index];
}//l
//				printf("\n%i %i  %i %i\n",d,dd,ddd,dddd);
				if(d ==3 && dd==2 && ddd==1 && dddd==0){
					  for(int i=0; i <len; i++){
         				       sixD[i] = pot6[i];
       						 }
				for(int c = 5; c >1; c--){
					for(int cc = c-1;cc>0;cc--){
						for(int ccc = cc-1;ccc>-1;ccc--){
							for(int n = 0; n <len; n++){
								int new_index = 0;
								for(int u =0;u<6;u++){
								stored_values[n][u] = dummy[n][u];}
								
								stored_values[n][c] = 0;
								stored_values[n][cc]=0;
								stored_values[n][ccc]=0;
								for(int k = 0; k<6;k++){
									new_index = new_index + stored_values[n][k];
								}
								new_index = new_index + corr;
								pot6[index[n]+corr] = pot6[index[n]+corr]-sixD[new_index];
}//n
//							printf("\n%i %i %i\n",c,cc,ccc);

				if(c==2 && cc==1 && ccc==0){
					  for(int i=0; i <len; i++){
                                               sixD[i] = pot6[i];
                                                 }

				for(int b = 5; b>0;b--){
					for(int bb = b-1;bb>-1;bb--){
						for(int p = 0; p<len;p++){
							int new_index = 0;
						for(int u = 0; u<6;u++){
							stored_values[p][u] = dummy[p][u];
						}
						stored_values[p][b] = 0;
						stored_values[p][bb] = 0;
						for(int k = 0; k<6;k++){
							new_index = new_index + stored_values[p][k];
						}
						new_index = new_index + corr;		
						pot6[index[p]+corr] = pot6[index[p]+corr]-sixD[new_index];
}//p
						//printf("\n%i %i\n",b,bb);
				if(b==1 && bb==0){
					  for(int i=0; i <len; i++){
                                               sixD[i] = pot6[i];
                                                 }
					for(int a = 5; a>-1;a--){
						for(int q = 0; q<len;q++){
							int new_index = 0;
						for(int u = 0; u<6; u++){
							stored_values[q][u] = dummy[q][u];
						}
						stored_values[q][a] = 0;
						for(int k = 0; k <6; k++){
							new_index = new_index + stored_values[q][k];
						}
						new_index = new_index + corr;
						pot6[index[q]+corr] = pot6[index[q]+corr]-sixD[new_index];
}//q

//						printf("\n%i\n",a);
}//a
}//if b,bb
}//bb
}//b
}//if c==2,cc,ccc
}//ccc
}//cc
}//c		
				
}//if d,dd,ddd,dddd
}//dddd
}//ddd
}//dd
}//d



}//if j==5
}//j


//for(int q = 0; q <len;q++){printf("\n%f %i\n",pot6[q],q);}

}//end method
//////////////////////////////////////////////////////////////////////////////


/*
//double cpot::get_coupling(double* pot,int degree,int grid_size,int dof){
void cpot::get_coupling(vector<double>& pot, int degree) {
if( degree ==6)
get6d(pot);
else if( degree ==5)
get5d(pot);
else if( degree == 4)
get4d(pot);
else if (degree == 3)
//get3d(pot,grid_size,dof);
get3d(pot);
else
//get2d(pot,grid_size,dof);
get2d(pot);
}

int cpot::fact(int n){
if(n==0 || n==1)
	return 1;
else
return n*fact(n-1);
}
*/
