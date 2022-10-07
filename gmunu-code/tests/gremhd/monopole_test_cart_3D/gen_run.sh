#!/bin/bash

#make allclean
make -j

LIMITER="woodward"

for KAPPA in "1" "5" "10" 
do
        echo "Generating KAPPA=${KAPPA}, LIMITER=${LIMITER} parameters file..."
	cp -r template_para.par k${KAPPA}_limiter_${LIMITER}.par
	cp -r mpi-condor.sub k${KAPPA}_limiter_${LIMITER}.sub
	sed -i "s/KAPPA/${KAPPA}/g"  k${KAPPA}_limiter_${LIMITER}.par
	sed -i "s/LIMITER/${LIMITER}/g" k${KAPPA}_limiter_${LIMITER}.par
	sed -i "s/para/k${KAPPA}_limiter_${LIMITER}/g" k${KAPPA}_limiter_${LIMITER}.sub
        condor_submit k${KAPPA}_limiter_${LIMITER}.sub
	
	#python get_norm.py
done
