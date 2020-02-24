#!/bin/bash
cd "${0%/*}"
module load python/3.8.0 scipy-stack
./compute_order_parameter.py $1