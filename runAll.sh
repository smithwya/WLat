#!/bin/bash
nJobs=1
beta=2.25
Nc=2
N=24
T=24
xiR=2
hot=1
nSweeps=0
gFixing=0
suffix=GRT
measurements=1
sweeps_between_meas=20
runtime=01:00:00
filepath=/N/project/Lattice-C/"SU("$Nc")"
#makes 'Configs' and 'Data' folders in filepath location



for beta in 2.75 #2.25 2.30 2.35 2.40 2.45 2.50 2.55 2.60 2.65 2.70 2.75
do
for xiR in 3 #1 2 3 4
do
make && bash runBatch.sh $nJobs $beta $N $(($xiR*$N)) $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $((6*$xiR)):00:00 $filepath $Nc
done
done

