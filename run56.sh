#!/bin/bash
nJobs=1
beta=6.8695638
Nc=3
N=56
T=56
xiR=1
hot=0
nSweeps=350
gFixing=0
suffix=GRT
measurements=0
sweeps_between_meas=20
runtime=24:00:00
filepath=/N/project/Lattice-C/SU"$Nc"
memory=18G

#makes 'Configs' and 'Data' folders in filepath location

make && bash runBatch.sh $nJobs $beta $N $(($xiR*$N)) $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $runtime $filepath $Nc $memory
