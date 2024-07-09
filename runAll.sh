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
measurements=10
sweeps_between_meas=200
runtime=24:00:00
filepath=/N/project/Lattice-C/"SU("$Nc")"
memory=4G

#makes 'Configs' and 'Data' folders in filepath location



for beta in 2.65 #2.50
do
for xiR in 8 #7 8
do
make && bash runBatch.sh $nJobs $beta $N $(($xiR*$N)) $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $((4*$xiR)):00:00 $filepath $Nc $(($xiR/2))G
done
done

