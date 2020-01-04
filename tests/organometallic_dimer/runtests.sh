#!/bin/bash

echo "********************"
echo "Different rho cutoff:"
echo "********************"
for dir in ./*/ ; do
   cd $dir
   echo "Time for execution in "${dir}
   time ../../../src_nciplot_4.0/nciplot < orgmet_dimer.nci > orgmet_dimer.nco
   int_n1=$(grep 'n=1.0' orgmet_dimer.nco | head -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   int_n1=$(echo ${int_n1} | bc -l ) 
   int_n1_r1=$(grep 'n=1.0' orgmet_dimer.nco | head -n 4 | tail -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   int_n1_r2=$(grep 'n=1.0' orgmet_dimer.nco | head -n 6 | tail -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   int_n1_r3=$(grep 'n=1.0' orgmet_dimer.nco | head -n 8 | tail -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   echo "Results from intramolecular calculation : "${int_n1}
   int_ranges=$(echo ${int_n1_r1} + ${int_n1_r2} + ${int_n1_r3} | bc -l )
   echo "                          and by ranges :" ${int_ranges}
   cd ..
done
rm -rf */*.cube
rm -rf */*.dat
rm -rf */*.vmd
rm -rf */*.nco
