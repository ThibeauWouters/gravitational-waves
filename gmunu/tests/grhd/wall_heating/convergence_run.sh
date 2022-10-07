#!/bin/bash

make allclean
make -j

for i in 32 64 128 256 512 1024 2048
do  
   echo "Generating NX=${i} parameters file..."
	cp -r template.par res_${i}.par
	sed -i "s/NX/${i}/g" res_${i}.par
        mpirun -np 8 ./gmunu -i res_${i}.par
done

