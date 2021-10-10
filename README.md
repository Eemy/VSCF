# VSCF
2-mode to 6-mode VSCF Solver

# Three Folders:
## CodeWriting ##
The key script is writeVSCF.csh here. It basically concatenates vscfFragment1.cpp, calls vscfWriter.cpp to write the body of the vscf code, then slaps on vscfFragment2.cpp at the end. It then compiles this final script called vscf.cpp and runs it. 

## Modules ##
Has code for the custom classes I wrote, as well as their header files for linking (Mode, EigSolver, Potential). cpot is the code written by Matt, which obtains differential coupling terms from Potentials. 

You will find that Potential.cpp is the chunkiest source file here. That's probably where most of the magic of my code takes place. Although it likely needs to be redesigned for better readability...

## UtilityFunc ##
Lapack libraries for diagonalization and matrix multiplication, other gauher routines, and AllocDouble are here. I didn't write these scripts, they're just inherited from other people that used them.

*IF YOU ARE LOOKING TO UNDERSTAND HOW ALL MY CODE COMES TOGETHER, YOU SHOULD READ THROUGH VSCF.CPP (sample file that was written and left outside of the folders) AND SEE HOW THE OBJECTS ARE USED IN THE MAIN BODY. THEN DIG INTO THE OBJECTS IN MY SOURCE CODE FOLDER AND FIGURE OUT WHAT THOSE FUNCTIONS ARE REALLY DOING UNDER THE HOOD. 

*THIS IS THE BEAUTY OF OOP, THERE IS A LOT OF ABSTRACTION IN THE CODE. VSCF.CPP SEEMS RELATIVELY SIMPLE AS IT SEEMS TO JUST BE INITIALIZING OBJECTS, LOOPING, AND CALLING FUNCTIONS. 


## IMPORTANT: Linking github to your shell ##

