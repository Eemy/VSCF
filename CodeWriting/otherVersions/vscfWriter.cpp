#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>

std::string makeLoopHeader(std::string var, std::string start, std::string end);
void print(FILE* script, std::string line);

int main(int argc, char** argv) {
  if(argc!=3) {
    printf("Error: <couplingDegree> <numMethods>\n");
  } else{
    int coupling = atoi(argv[1]);
    int numMethods = atoi(argv[2]);
    FILE *script = fopen("vscf.cpp","a");

    std::vector<int> dims; 
    if(numMethods > 1) {
      dims.resize(numMethods);
      std::cout << "Enter dim of potentials in ascending order\n"; 
      for(int i=0 ; i<numMethods ; i++) {
        std::cout << "Dim:";
        std::cin >> dims[i];
        if(dims[i] > coupling) {
          std::cout << "This dimension is not valid with the degree of coupling.\n";
          exit(0);
        }
      }
    }

    //Fill string array to hold loop variables
    std::string* loopVars = new std::string[coupling+1];
    for(int loopCount=0 ; loopCount<coupling+1 ; loopCount++) {
      char insert = 'a';
      insert += loopCount;
      std::string toString(1,insert);
      loopVars[loopCount] = insert;
    }
 
/////////////////////////////////////////////////VSCF.CPP/////////////////////////////////////////////////////////
    //Before iterations: Potential objects
    if(numMethods == 1) {
      print(script,"Potential pot(V,1,couplingDegree,nModes,nPoints,expectedLength,dof);\n");
      print(script,"Potential dx(Dx,1,couplingDegree,nModes,nPoints,expectedLength,dof);\n");
      print(script,"Potential dy(Dy,1,couplingDegree,nModes,nPoints,expectedLength,dof);\n");
      print(script,"Potential dz(Dz,1,couplingDegree,nModes,nPoints,expectedLength,dof);\n\n");
      print(script,"double** slices = pot.get1DSlices();\n\n");

      print(script,"std::vector<double**> dipSlices;\n");
      print(script,"dipSlices.push_back(dx.get1DSlices());\n");
      print(script,"dipSlices.push_back(dy.get1DSlices());\n");
      print(script,"dipSlices.push_back(dz.get1DSlices());\n");

      
      print(script,"std::vector<Potential*> dipoles;\n");
      print(script,"dipoles.push_back(&dx);\n");
      print(script,"dipoles.push_back(&dy);\n");
      print(script,"dipoles.push_back(&dz);\n");
    } else {
      for(int i=0 ; i<numMethods ; i++) {
        if(i==0) {
          std::string index = std::to_string(i);
          print(script,"Potential pot" + index + " (potentials[" + index + "],1,potDims[" + index + "],nModes,nPoints,expectedLengths[" + index + "],dof);\n");
          print(script,"double** slices = pot0.get1DSlices();\n");
        } else {
          std::string index = std::to_string(i);
          std::string index2 = std::to_string(i-1);
          print(script,"Potential pot" + index + " (potentials[" + index + "],potDims[" + index2 + "]+1,potDims[" + index + "],nModes,nPoints,expectedLengths[" + index + "],dof);\n");  
        }
      }   
    }
    print(script,"\n");

    print(script,"//Prepare: eigensolver on pure 1D slices for each mode\n");
    print(script,makeLoopHeader("i","0","nModes").c_str());
    print(script,"prevEnergy += solver.solveMode(dof[i],slices[i],0);\n}\n");
    
    //Ground state and Excited State Iterations        
    for(int i=0 ; i<2 ; i++) {
      if(i!=0) {
        fprintf(script,makeLoopHeader("z","0","nModes").c_str());
        print(script,"prevEnergy = 0.0;\n");
        fprintf(script,makeLoopHeader("i","0","nModes").c_str());
        print(script,"if(i==z) {\nprevEnergy += solver.solveMode(dof[i],slices[i],1);\n");
        print(script,"} else {\n prevEnergy += solver.solveMode(dof[i],slices[i],0);\n}\n}\n");
      }
      fprintf(script,makeLoopHeader("iter","1","maxIter").c_str());
      print(script,"int counter = 0;\n");  
      fprintf(script,makeLoopHeader("i","0","nModes").c_str());
      fprintf(script,makeLoopHeader("j","0","nPoints").c_str());
      print(script,"effV[i][j] = slices[i][j];\n}\n}\n");

//==================================AGH EFFECTIVE POTENTIALS===================================

      if(numMethods == 1) {
        //Loop through all the modes
        fprintf(script,makeLoopHeader(loopVars[0],"0","nModes").c_str());
        for(int loopCount = 1 ; loopCount < coupling ; loopCount++) {
          std::string start(loopVars[loopCount-1] + "+1");
          fprintf(script,makeLoopHeader(loopVars[loopCount],start,"nModes").c_str());
          std::ostringstream ss;
          ss << loopCount;
        }
      
        //Effective Potential Integral Loops
        fprintf(script,makeLoopHeader(loopVars[coupling],"0","nPoints").c_str());
        for(int q=0 ; q<coupling ; q++) {
          std::ostringstream ss;
          ss << q;
          print(script,("effV[" +loopVars[q]+"]["+loopVars[coupling]+"] += pot.integralDriver(counter, " + ss.str() + ", " + loopVars[coupling] +");\n"));
        }
        print(script,"}\ncounter++;\n");
        for(int loopCount = 0 ; loopCount < coupling ; loopCount++) {
          print(script,"}\n");
        }
      } else {
        for(int i=0 ; i<numMethods ; i++) {
          //Loop through all the modes
          print(script,"int counter" + std::to_string(i) + " = 0;\n");  
          fprintf(script,makeLoopHeader(loopVars[0],"0","nModes").c_str());
          for(int loopCount = 1 ; loopCount < dims[i] ; loopCount++) {
            std::string start(loopVars[loopCount-1] + "+1");
            fprintf(script,makeLoopHeader(loopVars[loopCount],start,"nModes").c_str());
            std::ostringstream ss;
            ss << loopCount;
          }
        
          //Effective Potential Integral Loops
          fprintf(script,makeLoopHeader(loopVars[dims[i]],"0","nPoints").c_str());
          for(int q=0 ; q<dims[i] ; q++) {
            std::ostringstream ss;
            ss << q;
            print(script,("effV[" +loopVars[q]+"]["+loopVars[dims[i]]+"] += pot" + std::to_string(i) + ".integralDriver(counter"+std::to_string(i)+", " + ss.str() + ", " + loopVars[dims[i]] +");\n"));
          }
          print(script,"}\ncounter"+std::to_string(i)+"++;\n");
          for(int loopCount = 0 ; loopCount < dims[i] ; loopCount++) {
            print(script,"}\n");
          }
        }
      }
//==================================================================================

      //Solve all 1D Schrodinger Eqns and sum energies
      int stateToGet = 0;
      if(i != 0) {
        stateToGet = i;
      }    
      print(script,"double energy = 0.0;\n");
      fprintf(script,makeLoopHeader("i","0","nModes").c_str());
      if(i != 0) {
        print(script,"if(i==z) {\n");
        std::ostringstream ss;
        ss << i;
        print(script,"energy += solver.solveMode(dof[i],effV[i]," + ss.str() + ");\n} else {\n"); 
      }
      print(script,"energy += solver.solveMode(dof[i],effV[i],0);\n");
      if(i != 0) {
        print(script,"}\n");
      }
      print(script,"}\n");

      //Get VSCF Energy
      if(numMethods == 1) {
        print(script,"energy -= pot.getVMinus();\n");
      } else {
        for(int i=0 ; i<numMethods ; i++) {
          print(script,"energy -= pot" + std::to_string(i) + ".getVMinus();\n");
        }
      }

      //Check for Convergence
      print(script,"if(checkConvergence(dof,energy,nModes)) {\n");
      print(script,"fprintf(results,\"Converged at iteration REPLACEd\\n\",iter);\n");
      if(i != 0) {
        print(script,"fprintf(results,\"Mode REPLACEi Excited-State VSCF Energy is: REPLACE.8f\\n\",z,energy*219474.6313708);\n");
        print(script,"excitedEnergies[z+1] = energy;\nbreak;\n} else {\n");
        print(script,"prevEnergy = energy;\n}\n");      
        print(script,"if(iter == maxIter-1) {\nfprintf(results,\"Mode REPLACEi VSCF failed to converge.\\n\",z);\nexcitedEnergies[z+1] = energy;\n}");
      } else {
        print(script,"fprintf(results,\"Ground-State VSCF Energy is: REPLACE.8f\\n\", energy*219474.6313708);\n");
        print(script,"excitedEnergies[0] = energy;\nbreak;\n} else {\n");
        print(script,"prevEnergy = energy;\n}\n");      
        print(script,"if(iter == maxIter-1) {\nprint(results,\"Ground-State VSCF failed to converge.\\n\");\nexcitedEnergies[0] = energy;\n}");

      }
     
      //Close remaining loops
      print(script,"\n}\n");
      if(i != 0) {
        print(script,"dof[z]->setExcitedState();\n}\n");
        print(script,"////////End Excited-State VSCF///////\n\n");
      } else {
        print(script,"///////End Ground-State VSCF/////////\n");
        fprintf(script,makeLoopHeader("i","0","nModes").c_str());
        print(script,"dof[i]->setGroundState();\n}\n");
        print(script,"/////////Excited-State VSCF//////////\n");
      }

    }//end i states

/////////////////////////////////////////DIPOLES////////////////////////////////////////////
    fprintf(script,makeLoopHeader("i","0","nModes").c_str());
    print(script,"dof[i]->updateWaveFcn(dof[i]->getGState());\noverlaps.push_back(dof[i]->getOverlapEG());\n}\n\n");
    print(script,"//DIPOLE CALCULATIONS\n");
    
    fprintf(script,makeLoopHeader("comp","0","3").c_str());
    print(script,"int counter = 0;\n");

    fprintf(script,makeLoopHeader("a","0","nModes").c_str());
    print(script,"intensityComponents[3*a+comp] += dipoles[comp]->integrateSlice(dof[a],dipSlices[comp][a],true);\n");//<1|D(1)|0> Integrals
    for(int loopCount = 1 ; loopCount < coupling ; loopCount++) {
      std::string start(loopVars[loopCount-1] + "+1");
      fprintf(script,makeLoopHeader(loopVars[loopCount],start,"nModes").c_str());
      std::ostringstream ss;
      ss << loopCount;

      if(loopCount == 1) {
        print(script,"intensityComponents[3*a+comp] += dipoles[comp]->integrateSlice(dof[b],dipSlices[comp][b],false)*overlaps[a];\n");        
        print(script,"intensityComponents[3*b+comp] += dipoles[comp]->integrateSlice(dof[a],dipSlices[comp][a],false)*overlaps[b];\n");        
      }
    }
    for(int i=0 ; i<coupling ; i++) {
      print(script,"intensityComponents[3*"+ loopVars[i] + "+comp] += dipoles[comp]->getDipole(counter,"+std::to_string(i)+");\n"); 
    }
    print(script,"counter++;");
    for(int i=0 ; i<coupling ; i++) {
      print(script,"\n}");
    }
    print(script,"\n}\n");

    fclose(script);
    
    delete[] loopVars;
  }//end else
}

std::string makeLoopHeader(std::string var, std::string start, std::string end) {
  return "for(int " + var + " = " + start + " ; " + var + "< " + end + " ; " + var + "++) {\n";
}

void print(FILE* script, std::string line) {
  fprintf(script,line.c_str());
}
