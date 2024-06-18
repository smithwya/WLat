#!/bin/bash
nJobs=1
beta=5.9
Nc=2
N=32
T=32
xiR=2
hot=0
nSweeps=3000
gFixing=0
suffix=GRT
measurements=0
sweeps_between_meas=0
runtime=24:00:00
filepath=/N/project/Lattice-C/"SU("$Nc")"
memory=4G

#makes 'Configs' and 'Data' folders in filepath location



for beta in 2.30
do
for xiR in 2
do
make && bash runBatch.sh $nJobs $beta $N $(($xiR*$N)) $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $((12*$xiR)):00:00 $filepath $Nc $((xiR*2))G
done
done

