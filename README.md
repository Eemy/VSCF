# VSCF
2-mode to 6-mode VSCF Solver

# Three Folders:
## CodeWriting ##
The key script is writeVSCF.csh here. It basically concatenates vscfFragment1.cpp, calls vscfWriter.cpp to write the body of the vscf code, then slaps on vscfFragment2.cpp at the end. It then compiles this final script called vscf.cpp and runs it. 
UPDATE: This part of the code was overhauled. Check the 'revamp' branch for more user-friendly code.

## Modules ##
Has code for the custom classes I wrote, as well as their header files for linking (Mode, EigSolver, Potential). cpot is the code written by Matt, which obtains differential coupling terms from Potentials. 

You will find that Potential.cpp is the chunkiest source file here. That's probably where most of the magic of my code takes place. Although it likely needs to be redesigned for better readability...

## UtilityFunc ##
Lapack libraries for diagonalization and matrix multiplication, other gauher routines, and AllocDouble are here. I didn't write these scripts, they're just inherited from other people that used them.

*IF YOU ARE LOOKING TO UNDERSTAND HOW ALL MY CODE COMES TOGETHER, YOU SHOULD READ THROUGH VSCF.CPP (sample file that was written and left outside of the folders) AND SEE HOW THE OBJECTS ARE USED IN THE MAIN BODY. THEN DIG INTO THE OBJECTS IN MY SOURCE CODE FOLDER AND FIGURE OUT WHAT THOSE FUNCTIONS ARE REALLY DOING UNDER THE HOOD. 

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
