#!/bin/bash

#make allclean
make -j

#number of cores
NP=4

#LIMITER="woodward"
#LIMITER="weno"

for TYPE_CT in "average" "uct_hll" "uct_contact" #"weno7" "exeno7" "teno5ad" "weno5cu6" 
do
	echo "Generating ct runs, TYPE_CT=${TYPE_CT} parameters file..."
	cp -r template_para.par ${TYPE_CT}.par
	sed -i "s/RUNNAME/${TYPE_CT}/g" ${TYPE_CT}.par
	sed -i "s/TYPE_CT/${TYPE_CT}/g" ${TYPE_CT}.par
	mpirun -np ${NP} ./gmunu -i ${TYPE_CT}.par
done
