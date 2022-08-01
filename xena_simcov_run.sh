#!/bin/bash
#SBATCH --partition=singleGPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=0-01:00:00
#SBATCH --job-name=simcov_test
#SBATCh --output=simcov.out

module load upcxx/2020.10.0-tsz6
module load cmake/3.19.2-dd55

upcxx-run -n $SLURM_NTASKS -N $SLURM_JOB_NUM_NODES -- install/bin/simcov --config=covid_default_small.config --output=results_DES
