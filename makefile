#!/bin/csh

#icpc vscfFragment.cpp Mode.cpp EigSolver.cpp ghQuad.cpp jacobi.cpp Potential.cpp matmultAB.cpp eigsrt.cpp diagon.cpp -o vscfFragment 
g++ -std=c++11 -g -ggdb3 vscf3d.cpp Mode.cpp EigSolver.cpp ghQuad.cpp jacobi.cpp Potential.cpp matmultAB.cpp eigsrt.cpp diagon.cpp cpot.cpp -o vscf3d
#g++ -std=c++11 -g -ggdb3 vscf.cpp Mode.cpp EigSolver.cpp ghQuad.cpp jacobi.cpp Potential2.cpp matmultAB.cpp eigsrt.cpp diagon.cpp cpot.cpp -o vscf 
