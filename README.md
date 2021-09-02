# VSCF
2-mode to 6-mode VSCF Solver

## Three Folders: ##
# Code Writing
The key script is writeVSCF.csh here. It basically concatenates vscfFragment1.cpp, calls vscfWriter.cpp to write the body of the vscf code, then slaps on vscfFragment2.cpp at the end. It then compiles this final script called vscf.cpp and runs it. 

# Source Code
Has code for the custom classes I wrote, as well as their header files for linking (Mode, EigSolver, Potential)

# Utility Functions
jacobi.cpp for diagonalizing (this should probably get replaced by another diag routine) and other gauher routines, as well as AllocDouble are here. I didn't write these scripts, they're just inherited from other people that used them.

*IF YOU ARE LOOKING TO UNDERSTAND HOW ALL MY CODE COMES TOGETHER, YOU SHOULD READ THROUGH VSCF.CPP AND SEE HOW THE OBJECTS ARE USED IN THE MAIN BODY. THEN DIG INTO THE OBJECTS IN MY SOURCE CODE FOLDER AND FIGURE OUT WHAT THOSE FUNCTIONS ARE REALLY DOING UNDER THE HOOD. 
*THIS IS THE BEAUTY OF OOP, THERE IS A LOT OF ABSTRACTION IN THE CODE. VSCF.CPP SEEMS RELATIVELY SIMPLE AS IT SEEMS TO JUST BE INITIALIZING OBJECTS, LOOPING, AND CALLING FUNCTIONS. 
