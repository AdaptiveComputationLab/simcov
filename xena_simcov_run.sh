#!/bin/bash
#SBATCH --partition=singleGPU
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --time=0-01:00:00
#SBATCH --job-name=simcov_test
#SBATCh --output=simcov.out

module load upcxx/2020.10.0-python3-tsz6
module load cmake/3.18.4-2lmi

upcxx-run -n $SLURM_NTASKS -N $SLURM_JOB_NUM_NODES -- install/bin/simcov --config=covid_default.config --output=results
