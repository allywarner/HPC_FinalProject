#!/bin/bash
#SBATCH -J mpi_partition
#SBATCH -o partition.o%j
#SBATCH -N 8 -n 8
#SBATCH -t 00:30:00
#SBATCH --mail-user=avw0731@gmail.com
#SBATCH --mail-type=end
#SBATCH --account=owner-guest
#SBATCH --partition=ember-guest

module load mpich2

mpirun ./partition Matrices/p_finance256.dat 100 1
