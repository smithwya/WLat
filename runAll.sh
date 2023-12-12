#!/bin/bash
nJobs=100
beta=6.0
Nc=3
N=24
T=96
xiR=5
hot=1
nSweeps=500
gFixing=1e-7
suffix=GRT
measurements=0
sweeps_between_meas=50
runtime=24:00:00
filepath=/N/slate/smithwya/"SU(3)"
#makes 'Configs' and 'Data' folders in filepath location

for xiR in 1 2 3 4
do
make && bash runBatch.sh $nJobs $beta $N $T $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $runtime $filepath $Nc
done
