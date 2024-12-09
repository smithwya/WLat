#!/bin/bash
Nc=3
beta=6.06293868
L=4
T=4
xiR=1

multigen=0
nTherm=0
gFixTol=1e-7

measure=1
Rmax=2
Tmax=2

#filepath=/N/project/Lattice-C/SU"$Nc"/Configs/"$beta"-"$xiR" 
#datpath=/N/project/Lattice-C/SU"$Nc"/Data/"$beta"-"$xiR"/GRT
startcnfg=200
endcnfg=400
inc=20

runtime=24:00:00
memory=2G


#mkdir -p $filepath
#mkdir -p $datpath


#make && sbatch submit.script --time=$runtime --mem=$memory $filepath/run $datpath/run $startcnfg $endcnfg $inc $Nc $beta $L $T $xiR $multigen $nTherm $gFixTol $measure $Rmax $Tmax

./bin/WLat run datrun $startcnfg $endcnfg $inc $Nc $beta $L $T $xiR $multigen $nTherm $gFixTol $measure $Rmax $Tmax
