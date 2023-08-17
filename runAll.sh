#!/bin/bash
nJobs=100
beta=2.6
N=20
T=20
xiR=1
hot=0
nSweeps=1000
gFixing=1e-7


bash runBatch.sh $nJobs $beta $N $T $xiR $hot $nSweeps $gFixing

