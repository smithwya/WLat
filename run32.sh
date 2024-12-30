#!/bin/bash
Nc=3
beta=6.4375
L=32
T=32
xiR=1

multigen=1
nTherm=0
gFixTol=0

measure=0
Rmax=10
Tmax=10

filepath=/N/project/Lattice-C/SU"$Nc"/Configs/"$L"^3x"$T"/"$beta"-"$xiR" 
datpath=/N/project/Lattice-C/SU"$Nc"/Data/"$L"^3x"$T"/"$beta"-"$xiR"/GRT
startcnfg=6720
endcnfg=7000
inc=20

runtime=04:00:00
memory=4G


mkdir -p $filepath
mkdir -p $datpath


#make && sbatch submit.script --time=$runtime --mem=$memory $filepath/run $datpath/run $startcnfg $endcnfg $inc $Nc $beta $L $T $xiR $multigen $nTherm $gFixTol $measure $Rmax $Tmax

sbatch --time=$runtime --mem=$memory submit.script $filepath/run $datpath/GRT"$i" $startcnfg $endcnfg $inc $Nc $beta $L $T $xiR $multigen $nTherm $gFixTol $measure $Rmax $Tmax
