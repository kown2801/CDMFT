#!/bin/bash
module reset
module load gcc nixpkgs boost
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/local/lib:${HOME}/local/lib/openblas
echo $LD_LIBRARY_PATH
