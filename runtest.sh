#!/bin/bash
nJobs=6
beta=6.06293868
Nc=3
N=28
T=28
xiR=1
hot=1
nSweeps=0
gFixing=1e-7
suffix=GRT
measurements=10
sweeps_between_meas=20
runtime=24:00:00
filepath=/N/project/Lattice-C/SU"$Nc"
memory=2G

#makes 'Configs' and 'Data' folders in filepath location

make && bash runBatch.sh $nJobs $beta $N $(($xiR*$N)) $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $runtime $filepath $Nc $memory
