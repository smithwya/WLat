#!/bin/bash
Nc=3
beta=6.2264790
L=24
T=24
xiR=1

multigen=1
nTherm=0
gFixTol=1e-7

measure=1
Rmax=12
Tmax=12

filepath=/N/project/Lattice-C/SU"$Nc"/Configs/"$L"^3x"$T"/"$beta"-"$xiR" 
datpath=/N/project/Lattice-C/SU"$Nc"/Data/"$L"^3x"$T"/"$beta"-"$xiR"/GRT
startcnfg=3000
endcnfg=7000
inc=20

runtime=24:00:00
memory=4G


mkdir -p $filepath
mkdir -p $datpath


#make && sbatch submit.script --time=$runtime --mem=$memory $filepath/run $datpath/run $startcnfg $endcnfg $inc $Nc $beta $L $T $xiR $multigen $nTherm $gFixTol $measure $Rmax $Tmax


#sbatch --time=$runtime --mem=$memory submit.script $filepath/run $datpath/GRT"$i" $i $(($i+80)) $inc $Nc $beta $L $T $xiR $multigen $nTherm $gFixTol $measure $Rmax $Tmax

for i in {3000..7000..100}
do
sbatch --time=$runtime --mem=$memory submit.script $filepath/run $datpath/GRT"$i" $startcnfg $endcnfg $inc $Nc $beta $L $T $xiR $multigen $nTherm $gFixTol $measure $Rmax $Tmax
done
