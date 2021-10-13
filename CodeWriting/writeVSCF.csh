#!/bin/csh

if($#argv != 4) then
  echo "<NModes> <NPoints> <Coupling> <NumMethods>"
  exit
endif

set origDir = `pwd`
#This line changes depending on where this code writing script is located
cd ~/git/VSCF/CodeWriting

cat vscfFragment1.cpp > vscf.cpp
vscfWriter $3 $4
sed -i "s/REPLACE/%/g" vscf.cpp

cat vscfFragment2.cpp >> vscf.cpp 

#RUN THE PROGRAM AND DELETE IT
echo "Compiling..."
make
cd $origDir
echo "Running vscf.cpp..."
#~/git/VSCF/CodeWriting/vscf $1 $2 $3 $4
#rm -f ~/git/VSCF/CodeWriting/vscf.cpp
#rm -f ~/git/VSCF/CodeWriting/vscf
echo "Done!"
