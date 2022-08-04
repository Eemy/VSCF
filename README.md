# VSCF
2-mode to 6-mode VSCF Solver

# Three Folders:
## CodeWriting ##
**UPDATE 8/3/2022:** The revamped version does not require any code writing. The only script needed is vscf.cpp with the executable vscf. It can be compiled with the makefile. Some major changes involving argument handling need to be understood to use this version of the code.

-Command line arguments lead with \<Nmodes\> and \<Npoints\>. Then after that point, you enter a series of 5 commands in this order: \<Name of Energy File\> \<Name of Dipole X File\> \<Name of Dipole Y File\> \<Name of Dipole Z File\> \<Dimensionality of Files\>. 

-The importance of this revamp is that you can enter any number of these set of 5 commands after the first one. For example, if the first set of files is for the 2MR, and you only want to paint in a couple of the triples beyond that, you may include just those triples in another set of files and enter them in another set of 5. This removes prep work needed to be done by the user.  

-Illustrate with example: Let's say your 2D scans give V2.dat, D2x.dat, D2y.dat, and D2z.dat. And then you decide to run scans on 10 triples of your choice. Those would go in V3.dat, D3x.dat, D3y.dat, D3z.dat. To call the code, you simply input: 
\<path\>/vscf \<nModes\> \<nPoints\> V2.dat D2x.dat D2y.dat D2z.dat 2 V3.dat D3x.dat D3y.dat D3z.dat 3 [...if more files, include args]

-What matters is that the 1st set of 5 arguments includes the scans you want to extract your 1D slices from. 

-LASTLY, more files need to be included upfront so the code knows which tuples you are representing in the files. From the previous example, you need files named "2.dat" and "3.dat" in the same directory of your potentials and freq.dat files. See sample files in this directory on how they should be formatted. 

Outdated: ~~The key script is writeVSCF.csh here. It basically concatenates vscfFragment1.cpp, calls vscfWriter.cpp to write the body of the vscf code, then slaps on vscfFragment2.cpp at the end. It then compiles this final script called vscf.cpp and runs it.~~ 

## Modules ##
Has code for the custom classes I wrote, as well as their header files for linking (Mode, EigSolver, Potential). cpot is the code written by Matt Laskowski, which obtains differential coupling terms (2d to 6d) from Potentials. 

You will find that Potential.cpp is the chunkiest source file here. That's probably where most of the action of my code takes place. Although it likely needs to be redesigned for better readability...

## UtilityFunc ##
Lapack libraries for diagonalization and matrix multiplication, other gauher routines, and AllocDouble are here. I didn't write these scripts, they're just inherited from other people that used them.

*IF YOU ARE LOOKING TO UNDERSTAND HOW ALL MY CODE COMES TOGETHER, YOU SHOULD READ THROUGH VSCF.CPP AND SEE HOW THE OBJECTS ARE USED IN THE MAIN BODY. THEN DIG INTO THE OBJECTS IN MY SOURCE CODE FOLDER AND FIGURE OUT WHAT THOSE FUNCTIONS ARE REALLY DOING UNDER THE HOOD. 

## IMPORTANT: Linking github to CHPC ##
You will need to add an SSH key to your github account
1. Type on cmd line: ssh-keygen -t rsa -b 4096 -C "\<UNID\>@kingspeak.chpc.utah.edu"
2. Press enter when asked where to save the key and enter a password
3. Type: eval \`ssh-agent -c\`
4. Type: ssh-add
5. Go to ~/.ssh/id_rsa.pub and copy the file's contents
6. Go to github >> Settings >> SSH and GPG Keys >> New SSH Key
7. Title it, then paste contents from .pub file into key box. Then click Create.
8. The next time you login, key permissions will be weird. You will need to change permissions by typing: \'chmod  400 ~/.ssh/id_rsa\' in order to push or pull between your local and this remote repository.

- Then to get the VSCF Repository on your shell, go to the repo page on web browser, click on the green button that says "Code", go to the "SSH" tab and copy the github URL.
- Go into your shell, type "git clone \<pasted URL from github\>"
- This should automatically upload all the folders and files from the repository.
- chmod +x all the files and you should be able to execute them! Learn some of the basic git commands so you can also modify the code and easily sync it up to the repository so everyone else can access the changes.
