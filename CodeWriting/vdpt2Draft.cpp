//======================VMP2 Corrections=======================
  std::vector<double> mp2Corr(nModes);

  //perturbation is total - effV
  int maxState = 3;
  int minState = 1;
  int numStates = maxState-minState+1;

///This way of writing singles misses certain integrals when the bra and ket of the excited mode are at the same quanta.
/*
  //Perturbation integrals involving single excitations 
  std::vector<double> singles(nModes*nModes*numStates);
  for(int i=0 ; i<potIterators.size() ; i++) {
    for(int j=0 ; j<potIterators[i].size() ; j++) {
      //Mode fixed at (1,0)
      for(int k=0 ; k<nModes ; k++) {
        dof[k]->setBra(1);
        //Excited Mode
        for(int l=0 ; l<potIterators[i][j].size() ; l++) {
          int secondMode = potIterators[i][j][l];
          //Quanta of the Single Excitation
          for(int m=minState ; m<=maxState ; m++) { //do not include |10..>
            dof[secondMode]->setKet(m); 
            double integralVal = pot[i]->integrateTuple(j,false);
            singles[firstMode*nModes*numStates+secondMode*numStates+(m-minState)] += integralVal;
          } 
          dof[secondMode]->setKet(0);
        }
        dof[k]->setBra(0);
      }
    }
  }
*/

//Another way to write singles
  std::vector<double> singles(nModes*nModes*numStates);
  for(int i=0 ; i<nModes ; i++) {
    dof[i]->setBra(1);
    for(int j=0 ; j<nModes ; j++) {
      std::vector<int> diff;
      diff.push_back(i);
      diff.push_back(j);
      for(int l=0 ; l<potIterators.size() ; l++) {
        for(int l2=0 ; l2<potIterators[l].size(); l2++) {
        
          for(int m=minState; m<=maxState ; m++) {
            dof[j]->setKet(m);

            double integralVal = 1.0;
            for(int n=0 ; n<diff.size() ; n++) {
              if(dof[diff[n]]->getBra() != dof[diff[n]]->getKet()) {
                bool found = false;
                for(int n2=0 ; n2<potIterators[l][l2].size() ; n2++) {
                  if(diff[n]==potIterators[l][l2][n2]) {
                    found = true;
                    break;
                  }
                }
                if(!found) {
                  integralVal = 0;
                  break;
                }
              }          
            }
            if(integralVal != 0) { 
              printf("Integral: %i %i %i, Tuple Num: %i\n",i,j,m,l2);        
              integralVal *= pot[l]->integrateTuple(l2,false);
              singles[i*nModes*numStates+j*numStates+(m-minState)] += integralVal;
            }
            dof[j]->setKet(0);
          }//m
        }//l2
      }//l
    }//j
    dof[i]->setBra(0);
  }//i
////////////////////////////////////////////////////////////////////

  //Perturbation integrals involving double excitations
  double nPairs = nModes*(nModes-1)/2;
  std::vector<double> doubles(nModes*nPairs*numStates*numStates);
  std::vector<int> indices(3); //hold tuple indices

  for(int i=0 ; i<nModes ; i++) {
    dof[i]->setBra(1);
    for(int j=0 ; j<nModes ; j++) {
      indices(0) = j;
      for(int j2=j+1 ; j2<nModes ; j2++) {
      indices(1) = j2;
        for(int j3=j2+1 ; j3<nModes ; j3++) {
          indices(2) = j3;

        for(int l=0 ; l<potIterators[i].size() ; l++) {
          for(int l2=0 ; l2<potIterators[i][j].size(); l2++) {
          
            for(int m=minState; m<=maxState ; m++) {
              dof[j]->setKet(m);
              for(int m2=minState ; m2<=maxState ; m2++) {
                dof[j2]->setKet(m2);
                for(int m3=minState ; m3<=maxState ; m3++) {              
                  dof[j3]->setKet(m3); 

                  double integralVal = 1.0;
                  for(int n=0 ; n<nModes ; n++) {
                    bool found = false;
                    for(int n2=0 ; n2<potIterators[i][j].size() ; n2++) P
                      if(n==potIterators[i][j][n2])
                        found = true;
                    }
                    if(!found)
                      integralVal *= dof[n]->getOverlapEG();
                    if(integralVal == 0)
                      break;
                  }

                  if(integralVal != 0) { 
                    double integralVal = pot[l]->integrateTuple(l2,false);
                    doubles[i*nPairs*numStates*numStates+
                            tupleIndexDriver(indices,nModes)*numStates*numStates+
                            (m-minState)*numStates+
                            (m2-minState)] += integralVal;
                  }

                  dof[j3]->setKet(0);
                }
                dof[j2]->setKet(0);
              }
              dof[j]->setKet(0);
            } 

          } 
        }
        
        }
      }
    }

  }





  for(int i=0 ; i<potIterators.size() ; i++) {
    for(int j=0 ; j<potIterators[i].size() ; j++) {
      //Mode fixed at (1,0)
      for(int k=0 ; k<potIterators[i][j].size() ; k++) {
        int firstMode = potIterators[i][j][k];
        dof[firstMode]->setBra(1);
        //Excited Mode
        for(int l=0 ; l<potIterators[i][j].size() ; l++) {
          //Excited Mode 2
          int secondMode = potIterators[i][j][l];
          indices[0] = secondMode; 
          for(int l2=l+1 ; l2<potIterators[i][j].size() ; l2++) {
            int thirdMode = potIterators[i][j][l2]; 
            indices[1] = thirdMode;
            //Quanta of the Excitation
            for(int m=minState ; m<=maxState ; m++) {
              dof[secondMode]->setKet(m); 
              //Quanta of 2nd Excitation
              for(int m2=minState ; m2<=maxState ; m2++) {
                dof[thirdMode]->setKet(m2);
                double integralVal = pot[i]->integrateTuple(j,false);
                doubles[firstMode*nPairs*numStates*numStates+
                        tupleIndexDriver(indices,nModes)*numStates*numStates+
                        (m-minState)*numStates+
                        (m2-minState)] += integralVal;
                dof[thirdMode]->setKet(0);
              }
              dof[secondMode]->setKet(0);
            } 
          }
        }
        dof[firstMode]->setBra(0);
      }
    }
  }

int stateDifference() {
  int counter = 0;
  for(int i=0 ; i<nModes ; i++) {
    if(dof[i]->getBra() != dof[i]->getKet())
      counter++
  }
  return counter;
}
