#!/bin/bash
#SBATCH --time=2-00:00:00
##SBATCH --ntasks=2
#SBATCH --nodes=4
##SBATCH --mem-per-cpu=1G
#SBATCH --mem=31G
#SBATCH --ntasks-per-node=24
#SBATCH --account=def-tremblay
#SBATCH --job-name="dmft_supra"

module reset
module load gcc nixpkgs boost
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/local/lib:${HOME}/local/lib/openblas

./launch.py 20
