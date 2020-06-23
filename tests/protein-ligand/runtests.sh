#!/bin/bash

echo "********************"
echo "Checking accuracy:"
echo "********************"
for dir in ./*/ ; do
   cd $dir
   echo "Time for execution with the accuracy level of "${dir}
   time ../../../src_nciplot_4.0/nciplot < 2v5x.nci > 2v5x.nco
   int_n1=$(grep 'n=1.0' 2v5x.nco | head -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   int_n1=$(echo ${int_n1} | bc -l ) 
   echo "Results from integration in "${dir}" : "${int_n1}
   cd ..
done
rm -rf */*.cube
rm -rf */*.dat
rm -rf */*.vmd
rm -rf */*.nco
