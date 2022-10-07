#!/bin/bash

#make allclean
make -j

#number of cores
NP=8

#LIMITER="woodward"
#LIMITER="weno"

for LIMITER in "weno5" # "wenoyc3" #"weno5" "weno7" "exeno7" "teno5ad" "weno5cu6" 
do
	rm output_*
	rm res_*
	
	for i in 32 64 128 256 512 1024 2048
	do  
           j=$(( i/NP ))
	   echo "Generating NXB=${j}, NX=${i}, LIMITER=${LIMITER} parameters file..."
		cp -r template_para.par res_${i}_limiter_${LIMITER}.par
		sed -i "s/NXB/${j}/g" res_${i}_limiter_${LIMITER}.par
		sed -i "s/NX/${i}/g" res_${i}_limiter_${LIMITER}.par
		sed -i "s/LIMITER/${LIMITER}/g" res_${i}_limiter_${LIMITER}.par
		mpirun -np ${NP} ./gmunu -i res_${i}_limiter_${LIMITER}.par
	done
	
	python get_norm.py
done
