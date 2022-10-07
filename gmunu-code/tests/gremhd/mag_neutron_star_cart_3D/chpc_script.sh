#!/bin/bash
#SBATCH -J BU0_2D
#SBATCH -N 1 -c 32
#SBATCH --output=res.txt

module load intel/2020u1
#/opt/share/openmpi-3.0.1/bin/mpirun ./gmunu -i gmunu.par
srun ./gmunu -i gmunu.par
