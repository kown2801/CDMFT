#!/bin/bash
rm slurm*
rm logfile
cd IN/
rm *
cd ../OUT/
rm *
cd ../DATA/
rm *
cd ../
cd BACKUP_START
cp params0.meas.json ../OUT/
cp LinkN.json ../IN/
cp LinkA.json ../IN/
cp self.dat ../DATA/
