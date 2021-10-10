#!/bin/csh

if($#argv != 4) then
  echo "<NModes> <NPoints> <Coupling> <NumMethods>"
  exit
endif

set origDir = `pwd`
#This line changes depending on where this code writing script is located
cd ~/anharmCalc/eigensolver/personalSolvers/nModeVSCF/original9.30

cat vscfFragment1.cpp > vscf.cpp
vscfWriter $3 $4
sed -i "s/REPLACE/%/g" vscf.cpp

cat vscfFragment2.cpp >> vscf.cpp 

#RUN THE PROGRAM AND DELETE IT
echo "Compiling..."
#makefile
make
cd $origDir
echo "Running vscf.cpp..."
#~/anharmCalc/eigensolver/personalSolvers/nModeVSCF/vscf $1 $2 $3 $4
#rm -f ~/anharmCalc/eigensolver/personalSolvers/nModeVSCF/vscf.cpp
#rm -f ~/anharmCalc/eigensolver/personalSolvers/nModeVSCF/vscf
echo "Done!"
