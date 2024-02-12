#!/bin/bash
nJobs=100
beta=5.9
Nc=3
N=24
T=24
xiR=3
hot=1
nSweeps=0
gFixing=1e-7
suffix=GRT
measurements=5
sweeps_between_meas=50
runtime=24:00:00
filepath=/N/project/Lattice-C/"SU(3)"
#makes 'Configs' and 'Data' folders in filepath location

for beta in 5.9 6.1 6.2 6.3
do
make && bash runBatch.sh $nJobs $beta $N $(($xiR*$N)) $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $runtime $filepath $Nc

done
