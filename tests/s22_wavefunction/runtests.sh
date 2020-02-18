#!/bin/bash

# Declare some array variables
declare -a arr_ref=()
declare -a dirlist=()
declare -a intvals=()

echo "********************"
echo "Checking S22 dataset"
echo "********************"
# Run!
rm -rf correlation.dat
n=0
cat reference_test.txt |
while read line ; do
   case $line in 
   	   (*"wfx"*)
   	   id=${line%%wfx*}
	   gid=${line%%.wfx*}
   	   name=$gid".wfx"
   	   ref=${line##*wfx}
   	   arr=$name
	   arr_a="A"$name
	   arr_b="B"$name
	   linemerge=$(echo $line | tr -s " " )
	   energy=$(echo $linemerge | cut -d ' ' -f 2 )
	   
	   sed -ie '2s/.*/'$arr_a'/' test.nci
	   sed -ie '3s/.*/'$arr_b'/' test.nci

	   mkdir ${gid}
	   cp test.nci ${gid}/${gid}.nci
	   cp plot_2d.gp ${gid}/
	   cp dimers/${name} ${gid}/${arr}
	   cp monomers_A/${name} ${gid}/${arr_a}
	   cp monomers_B/${name} ${gid}/${arr_b}
	   cd ${gid}
	   sed -i -- "s/PLACEHOLDER/A${gid}/g" *.gp 
           echo "Time for execution of "${gid}
	   time ../../../src_nciplot_4.0/nciplot < ${gid}.nci > ${gid}.nco
	   int=$(grep 'n=2.0' ${gid}.nco | head -n 1 | tail -n 1 |  tr -s " " | cut -d ' ' -f 4 ) #n=2.0 works best
           gnuplot plot_2d.gp
	   cd ..
	   n=$(($n+1))
           intvals[$n]=${int}
	   arr_ref[$n]=${ref}
	   dirlist[$n]=${gid}
	   echo ${dirlist[$n]} ${arr_ref[$n]} ${intvals[$n]} >> correlation.dat
	   ;; 
   esac
done
newlist=$( ls -d */)
echo "Directories created:"
echo ${newlist[@]}
gnuplot correlation.gp
echo "********************"
echo "Removing S22 files!!"
echo "********************"
# Cleanup
cat reference_test.txt |
while read line ; do
   case $line in 
   	   (*"wfx"*)
   	   id=${line%%wfx*}
	   gid=${line%%.wfx*}
   	   name=$gid".wfx"
   	   ref=${line##*wfx}
	   arr_a="A"$name
	   arr_b="B"$name
           cd ${gid}
           cp *.png ../plots/
           cd ..
	   rm -rf ${gid}
	   ;; 
   esac
done
rm -rf test.ncie 
