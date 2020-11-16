#!/bin/bash
#SBATCH --time=0-00:10:00
#SBATCH --job-name="actions"
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-tremblay-ac
#SBATCH --error=occupation-log.out
#SBATCH --output=occupation-log.out

module reset
module load gcc/7.3.0 nixpkgs boost scipy-stack openmpi/3.1.2
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/local/lib:${HOME}/local/lib/openblas
 #Set the current working directory to the directory this file is in
./actions.py $@