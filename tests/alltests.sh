#!/bin/bash

echo "**********************"
echo "Will run every test!!!"
echo "(This will take a bit)"
echo "**********************"
for dir in ./*/ ; do
   cd $dir
   echo "Executing tests in "${dir}
   ./runtests.sh 
   cd ..
done
echo "**********************"
echo "Hope everything works!"
echo "**********************"
