#!/bin/bash
nJobs=100
beta=6.0
Nc=3
N=24
T=96
xiR=1
hot=0
nSweeps=500
gFixing=0
suffix=wilson_loops
measurements=0
sweeps_between_meas=50
runtime=24:00:00
filepath=/N/slate/smithwya/"SU(3)"
#makes 'Configs' and 'Data' folders in filepath location


make && bash runBatch.sh $nJobs $beta $N $T $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $runtime $filepath $Nc

