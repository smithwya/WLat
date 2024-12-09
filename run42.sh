#!/bin/bash
nJobs=2
beta=6.33676665
Nc=3
N=42
T=42
xiR=1
hot=1
nSweeps=0
gFixing=1e-7
suffix=GRT
measurements=2
sweeps_between_meas=20
runtime=48:00:00
filepath=/N/project/Lattice-C/SU"$Nc"
memory=6G

#makes 'Configs' and 'Data' folders in filepath location

make && bash runBatch.sh $nJobs $beta $N $(($xiR*$N)) $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $runtime $filepath $Nc $memory
