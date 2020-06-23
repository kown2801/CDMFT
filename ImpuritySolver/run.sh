#!/bin/bash
#SBATCH --time=0-01:00:00
#SBATCH --job-name="dmft_test"
#SBATCH --ntasks=96
#SBATCH --nodes=3
#SBATCH --mem-per-cpu=1G
#SBATCH --account=def-tremblay

module reset
module load gcc/7.3.0 nixpkgs boost scipy-stack openmpi/3.1.2
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/local/lib:${HOME}/local/lib/openblas

srun IS IN/ OUT/ params1