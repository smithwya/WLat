#!/bin/bash
nJobs=100
beta=2.25
N=8
T=8
xiR=1
hot=0
nSweeps=200
gFixing=1e-7


bash runBatch.sh $nJobs $beta $N $T $xiR $hot $nSweeps $gFixing

