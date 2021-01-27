#!/bin/bash
#SBATCH --partition=dualGPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=2-00:00:00
#SBATCH --job-name=simcov_lung_model_test_01
#SBATCh --output=simcov.out

module load upcxx/2020.10.0-python3-3o75
module load cmake/3.18.4-2lmi

cd lungmodel
./lungmodel true
cd ..
cp lungmodel/*.dat .
upcxx-run -n 1 -N 1 -- install/bin/simcov --dim [300,300,300] --timesteps 1 --progress -v

