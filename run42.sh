#!/bin/bash
nJobs=1
beta=6.33676665
Nc=3
N=42
T=42
xiR=1
hot=0
nSweeps=1500
gFixing=0
suffix=GRT
measurements=0
sweeps_between_meas=20
runtime=24:00:00
filepath=/N/project/Lattice-C/SU"$Nc"
memory=16G

#makes 'Configs' and 'Data' folders in filepath location

make && bash runBatch.sh $nJobs $beta $N $(($xiR*$N)) $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $runtime $filepath $Nc $memory
