#!/bin/bash
nJobs=1
beta=2.5
Nc=2
N=8
T=8
xiR=2
hot=1
nSweeps=100
gFixing=1e-7
suffix=debug
measurements=1
sweeps_between_meas=50
runtime=2:00:00
filepath=/N/slate/smithwya/"SU(2)"
#makes 'Configs' and 'Data' folders in filepath location

for xiR in 2 
do
make && bash runBatch.sh $nJobs $beta $N $T $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $runtime $filepath $Nc
done
