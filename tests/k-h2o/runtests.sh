#!/bin/bash

echo "********************"
echo "Checking accuracy:"
echo "********************"
echo "Time for execution:"
time ../../src_nciplot_4.0/nciplot < kh2o6.nci > kh2o6.nco
int_n1=$(grep 'n=1.0' kh2o6.nco | head -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
int_n1=$(echo ${int_n1} | bc -l ) 
int_n1_r1=$(grep 'n=1.0' kh2o6.nco | head -n 4 | tail -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
int_n1_r2=$(grep 'n=1.0' kh2o6.nco | head -n 6 | tail -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
int_n1_r3=$(grep 'n=1.0' kh2o6.nco | head -n 8 | tail -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
echo "Results from intramolecular integration  : "${int_n1}
int_ranges=$(echo ${int_n1_r1} + ${int_n1_r2} + ${int_n1_r3} | bc -l )
echo "                           and by ranges :" ${int_ranges}
gnuplot kh2o6.gp
rm -rf *.cube
rm -rf *.dat
rm -rf *.vmd
rm -rf *.nco
