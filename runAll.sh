#!/bin/bash
nJobs=100
beta=2.75
N=8
T=8
xiR=4
hot=0
nSweeps=1000
gFixing=1e-7
suffix=green
measurements=10
sweeps_between_meas=100
runtime=00:20:00
filepath=/N/slate/smithwya/"SU(2)"
#makes 'Configs' and 'Data' folders in filepath location

for xiR in 1 2 3 4
do
bash runBatch.sh $nJobs $beta $N $(($xiR*$T)) $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $runtime $filepath
done
