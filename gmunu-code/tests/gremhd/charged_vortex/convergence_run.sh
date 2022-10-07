#!/bin/bash

#make allclean
make -j

#number of cores
NP=16

#LIMITER="woodward"
#LIMITER="weno"

rm output_*
rm res_*

for LIMITER in "woodward" "ppm" #"wenozp5" # "wenoyc3" #"weno5" "weno7" "exeno7" "teno5ad" "weno5cu6" 
do
	
	for i in 32 64 128 256 512 1024 2048 4096
	do  
           j=$(( 16 ))
	   echo "Generating NXB=${j}, NX=${i}, LIMITER=${LIMITER} parameters file..."
		cp -r template_para.par res_${i}_limiter_${LIMITER}.par
		sed -i "s/NXB/${j}/g" res_${i}_limiter_${LIMITER}.par
		sed -i "s/NYB/${j}/g" res_${i}_limiter_${LIMITER}.par
		sed -i "s/NX/${i}/g" res_${i}_limiter_${LIMITER}.par
		sed -i "s/NY/${i}/g" res_${i}_limiter_${LIMITER}.par
		sed -i "s/LIMITER/${LIMITER}/g" res_${i}_limiter_${LIMITER}.par

		#mpirun -np ${NP} ./gmunu -i res_${i}_limiter_${LIMITER}.par

		cp -r mpi-condor.sub res_${i}_limiter_${LIMITER}.sub
		sed -i "s/NP/${NP}/g" res_${i}_limiter_${LIMITER}.sub
		sed -i "s/PARA/"res_${i}_limiter_${LIMITER}.par"/g" res_${i}_limiter_${LIMITER}.sub
		condor_submit -batch-name ${i}_${LIMITER} res_${i}_limiter_${LIMITER}.sub
            wait
	done
        wait
	#python get_norm.py
done
