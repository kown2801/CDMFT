#!/bin/bash
module reset
module load StdEnv/2020 scipy-stack python
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/local/lib
