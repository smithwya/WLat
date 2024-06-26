#!/bin/bash
nJobs=100
beta=2.25
Nc=2
N=24
T=24
xiR=2
hot=1
nSweeps=0
gFixing=1e-7
suffix=GRT
measurements=3
sweeps_between_meas=20
runtime=24:00:00
filepath=/N/project/Lattice-C/"SU("$Nc")"
memory=4G

#makes 'Configs' and 'Data' folders in filepath location



for beta in 2.70 2.75
do
for xiR in 3 4
do
make && bash runBatch.sh $nJobs $beta $N $(($xiR*$N)) $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $((6*$xiR)):00:00 $filepath $Nc $((xiR*1))G
done
done

