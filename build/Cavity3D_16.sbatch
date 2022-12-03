#!/bin/bash
#SBATCH -J Cavity3D_17
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=mpp3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mail-type=end
#SBATCH --mail-user=vikas.kurapati@tum.de
#SBATCH --export=NONE

module load slurm_setup
module unload intel-mpi/2019-intel
module unload intel-mkl/2020
module unload intel-oneapi-compilers/2021.4.0
module load gcc/11.2.0
module load openmpi/4.1.2-gcc11
module load cmake
module load petsc/3.17.2-gcc11-ompi-real

mpiexec -n $SLURM_NTASKS ./NS-EOF-Runner ../ExampleCases/Cavity3D_17.xml