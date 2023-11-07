#!/bin/bash
nJobs=1
beta=6.0
Nc=3
N=12
T=12
xiR=1
hot=0
nSweeps=30
gFixing=0
suffix=green
measurements=0
sweeps_between_meas=100
runtime=03:00:00
filepath=/N/slate/smithwya/"SU(3)"
#makes 'Configs' and 'Data' folders in filepath location


make && bash runBatch.sh $nJobs $beta $N $T $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $runtime $filepath $Nc

