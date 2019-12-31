#!/bin/bash

echo "********************"
echo "Check default cases:"
echo "********************"
for dir in ./*defaults*/ ; do
   cd $dir
   ../../../src_nciplot_4.0/nciplot < AT.nci > AT.nco
   int_n1=$(grep 'n=1.0' AT.nco | head -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   int_n1=$(echo ${int_n1} | bc -l ) 
   echo "Results from integration in "${dir}" : "${int_n1}
   cd ..
done

echo "********************"
echo "Intermolecular case:"
echo "********************"
for dir in ./*inter/ ; do
   cd $dir
   echo "Time for execution with the accuracy level of "${dir}
   time ../../../src_nciplot_4.0/nciplot < AT.nci > AT.nco
   int_n1=$(grep 'n=1.0' AT.nco | head -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   int_n1=$(echo ${int_n1} | bc -l ) 
   echo "Results from intermolecular integration in "${dir}" : "${int_n1}
   cd ..
done

for dir in ./*inter_range/ ; do
   cd $dir
   ../../../src_nciplot_4.0/nciplot < AT.nci > AT.nco
   int_n1=$(grep 'n=1.0' AT.nco | head -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   int_n1=$(echo ${int_n1} | bc -l ) 
   int_n1_r1=$(grep 'n=1.0' AT.nco | head -n 4 | tail -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   int_n1_r2=$(grep 'n=1.0' AT.nco | head -n 6 | tail -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   int_n1_r3=$(grep 'n=1.0' AT.nco | head -n 8 | tail -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   echo "Results from intermolecular integration in "${dir}" : "${int_n1}
   int_ranges=$(echo ${int_n1_r1} + ${int_n1_r2} + ${int_n1_r3} | bc -l )
   echo "                                      and by ranges :" ${int_ranges}
   cd ..
done

echo "********************"
echo "Intramolecular case:"
echo "********************"
for dir in ./*intra/ ; do
   cd $dir
   echo "Time for execution with the accuracy level of "${dir}
   time ../../../src_nciplot_4.0/nciplot < AT.nci > AT.nco
   int_n1=$(grep 'n=1.0' AT.nco | head -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   echo "Results from intramolecular integration in "${dir}" : "${int_n1}
   cd ..
done

for dir in ./*intra_range/ ; do
   cd $dir
   ../../../src_nciplot_4.0/nciplot < AT.nci > AT.nco
   int_n1=$(grep 'n=1.0' AT.nco | head -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   int_n1=$(echo ${int_n1} | bc -l ) 
   int_n1_r1=$(grep 'n=1.0' AT.nco | head -n 4 | tail -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   int_n1_r2=$(grep 'n=1.0' AT.nco | head -n 6 | tail -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   int_n1_r3=$(grep 'n=1.0' AT.nco | head -n 8 | tail -n 1 | tr -s ' ' | cut -d ' ' -f 4 )
   echo "Results from intramolecular integration in "${dir}" : "${int_n1}
   int_ranges=$(echo ${int_n1_r1} + ${int_n1_r2} + ${int_n1_r3} | bc -l )
   echo "                                      and by ranges :" ${int_ranges}
   cd ..
done

rm -rf */*.cube
rm -rf */*.dat
rm -rf */*.vmd
rm -rf */*.nco
