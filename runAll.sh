#!/bin/bash
nJobs=100
beta=5.9
Nc=3
N=32
T=32
xiR=4
hot=1
nSweeps=1000
gFixing=1e-7
suffix=GRT
measurements=1
sweeps_between_meas=20
runtime=24:00:00
filepath=/N/project/Lattice-C/"SU("$Nc")"
#makes 'Configs' and 'Data' folders in filepath location

for beta in 5.9 6.0 6.1 6.2
do
make && bash runBatch.sh $nJobs $beta $N $(($xiR*$N)) $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $runtime $filepath $Nc

done

