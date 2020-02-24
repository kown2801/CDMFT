#SBATCH --ntasks=96
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-tremblay-ac

module reset
module load gcc nixpkgs boost scipy-stack
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/local/lib:${HOME}/local/lib/openblas

