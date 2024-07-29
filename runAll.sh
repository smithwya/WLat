#!/bin/bash
nJobs=100
beta=5.9
Nc=3
N=32
T=32
xiR=1
hot=1
nSweeps=1000
gFixing=0
suffix=GRT
measurements=0
sweeps_between_meas=0
runtime=24:00:00
filepath=/N/project/Lattice-C/"SU("$Nc")"
memory=4G

#makes 'Configs' and 'Data' folders in filepath location



for beta in 5.9 6.0
do
for xiR in 3 4
do
make && bash runBatch.sh $nJobs $beta $N $(($xiR*$N)) $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $((8*$xiR)):00:00 $filepath $Nc $(($xiR/2))G
done
done

