#!/bin/bash
nJobs=100
beta=2.25
N=8
T=8
xiR=4
hot=0
nSweeps=1000
gFixing=1e-7
suffix=test2_green
measurements=10
sweeps_between_meas=100
runtime=04:00:00


bash runBatch.sh $nJobs $beta $N $(($xiR*$T)) $xiR $hot $nSweeps $gFixing $suffix $measurements $sweeps_between_meas $runtime

