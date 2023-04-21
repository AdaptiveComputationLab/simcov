#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --time=2-00:00
#SBATCH --job-name=simcov-tutorial
#SBATCH --mail-user=
#SBATCH --mail-type=END

# Setup environment
export GASNET_MAX_SEGSIZE=128MB/P
export UPCXX_SHARED_HEAP_SIZE='128 MB'
export UPCXX_THREADMODE=seq
export UPCXX_CODEMODE=opt
export UPCXX_NETWORK=ibv

# Load modules
module purge
module load gcc/12.1.0-crtl
module load cmake/3.11.4-qkyj
module load openmpi/4.1.3-j6zb
module load upcxx/2020.10.0-6eh2

srun --mpi=pmi2 install/bin/simcov --config covid_default.config
