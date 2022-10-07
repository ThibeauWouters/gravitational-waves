#!/bin/bash
#SBATCH -J test 
#SBATCH --output=outfile_%j.txt
#SBATCH -e err_%j.err       # output error
#SBATCH -N 1 -n 32 --mem=MaxMemPerNode

module load intel/2020u4
export LD_LIBRARY_PATH=/scratch/s1/TjonnieLi/gmunu_libs/2020u4/mpiifort/hdf5/lib/:$LD_LIBRARY_PATH
export I_MPI_DEBUG=4
#ulimit -c unlimited
ulimit -a
srun --mpi=pmi2 gmunu -i full.par
