#!/bin/bash

echo "********************"
echo "Check default cases:"
echo "********************"
for dir in ./*/ ; do
   cd $dir
   echo "Time for execution with the accuracy level of "${dir}
   time ../../../src_nciplot_4.0/nciplot < PhenolDimer.nci > PhenolDimer.nco
   int_n1=$(grep 'n=1.0' PhenolDimer.nco | head -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   int_n1=$(echo ${int_n1} | bc -l ) 
   echo "Results from intramolecular integration in "${dir}" : "${int_n1}
   cd ..
done

rm -rf */*.dat
rm -rf */*.cube
rm -rf */*.nco
