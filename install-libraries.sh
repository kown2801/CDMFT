#!/bin/bash
#First create the 'local' folder structure
cd
mkdir local
cd local
mkdir lib include
#Then install OpenBLAS
git clone https://github.com/xianyi/OpenBLAS.git
cd OpenBLAS
git checkout v0.3.9
make -j4
cp libopenblas.a ~/local/lib
#Then install nlohmann/json
cd ~/local
git clone https://github.com/nlohmann/json.git
cp json/single_include/nlohmann/ include/nlohmann -r
