#!/bin/zsh

cnfg_in=testLat
cnfg_out=testLat
dat=testDat

Nc=3
beta=6.0
L=4
T=4
xiR=1
nSweeps=10
multiGen=1

cnfgStart=0
cnfgEnd=50
dSweeps=10

gFix=1
gTol=1e-7

measure=1
Rmax=2
Tmax=2

echo "Hello?"
./bin/WLat "$cnfg_in $cnfg_out $dat $Nc $beta $L $T $xiR $nSweeps $multiGen $cnfgStart $cnfgEnd $dSweeps $gFix $gTol $measure $Rmax $Tmax"
