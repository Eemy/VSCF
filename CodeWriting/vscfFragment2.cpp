  //Print out all the transition frequencies
  print(results,"************************************************************************\n");
  print(results," VSCF ground-state energy: (cm^-1) \n");
  fprintf(results," % -15.4f \n", excitedEnergies[0]*(219474.6313708));
  print(results," \n");
  print(results," Transitions: (cm^-1) \n");
  print(results," \n");
  print(results,"  Harmonic        VSCF            Intensity(km/mol)\n");
  for(int i=0; i<nModes ; i++) {
  fprintf(results," % -15.4f % -15.4f \n", freq[i],(excitedEnergies[i+1]-excitedEnergies[0])*(219474.6313708));
  }

  //DEALLOCATE
    delete[] freq;
    delete[] mass;
    for(int i=0 ; i<nModes ; i++) {
      delete[] effV[i];
      delete[] slices[i];
      delete dof[i];
    }
    delete[] effV;
    delete[] slices;
    delete[] dof;
    delete[] excitedEnergies;
  }
  return 0;
}

//=====================================OTHER METHODS==========================================
void readin(double* freq, double* mass, int N) {
  //read in frequencies 
  std::ifstream in("freq.dat",std::ios::in);
  if(!in) {
    printf("Error: freq.dat could not be opened\n");
    exit(0);
  }

  for(int i=0 ; i<N ; i++) in >> freq[i];
  in.close();
  
  //read in reduced masses
  in.open("rmass.dat");
  if(!in) {
    printf("Error: rmass.dat could not be opened\n");
    exit(0);
  }
  for(int i=0 ; i<N ; i++) in >> mass[i];
  in.close();
} 


void readPot(std::string fileName, double* potential, int length, int potLength) {
  //read in potential and subtract out equilibrium energy
  std::ifstream in(fileName,std::ios::in);
  if(!in) {
    printf("Error: %s could not be opened\n",fileName.c_str());
    exit(0);
  }
  double value;
  int index = 0;
  while(in >> value) {
    if(index < potLength) {
      potential[index++] = value;
    } else {
      printf("%s is too long.\n",fileName.c_str());
      exit(0);
    } 
  }
  if(index < potLength) {
    printf("%s is too short.\n",fileName.c_str());
    exit(0);
  }
  in.close();

  double eqEnergy = potential[(length-1)/2];
  for(int i=0 ; i<potLength ; i++) {
    potential[i] -= eqEnergy;
  }
} 

bool checkConvergence(Mode** dof, double energy, int nModes) {
  double diff = 0.0;
  for(int i=0 ; i<nModes ; i++) {
    double temp = dof[i]->computeMaxDiff();
    if(temp > diff)
      diff = temp;
  }
  return (diff < 1.0E-5) && (fabs(energy-prevEnergy)*219474.6313708 <0.5);
}

int fact(int n) {
  if(n==0 || n==1)
    return 1;
  else
    return n*fact(n-1);
}

void print(FILE* script, std::string line) {
  fprintf(script,line.c_str());
}

int lengthCheck(int nMode, int coupling) {
  int num = nMode;
  int couplingFac = fact(coupling);
  for(int i=1; i<coupling ; i++)
    num *= (nMode-i);
  return num/couplingFac;
}
