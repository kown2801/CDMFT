#!/bin/bash
#SBATCH --time=0-01:00:00
#SBATCH --job-name="dmft_test"
#SBATCH --ntasks=4
##SBATCH --nodes=3
#SBATCH --mem-per-cpu=1G
#SBATCH --account=def-tremblay


source ../scripts/export.sh ; srun IS IN/ OUT/ params1
