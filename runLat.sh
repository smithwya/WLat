#!/bin/bash
Nc=3
beta=6.612035
L=40
T=40
xiR=1

multigen=0
nTherm=0
gFixTol=1e-7
measure=1
Rmax=20
Tmax=20

startcnfg=3000
endcnfg=11000
inc=20

meas_per_job=2
runtime=72:00:00
memory=5G

filepath=/N/project/Lattice-C/SU"$Nc"/Configs/"$L"^3x"$T"/"$beta"-"$xiR"
datpath=/N/project/Lattice-C/SU"$Nc"/Data/"$L"^3x"$T"/"$beta"-"$xiR"/GRT

mkdir -p $filepath
mkdir -p $datpath


if [ "$multigen" -gt "0" ]; then

make && sbatch -o "Logs/L$L-b$beta-s%j-g$startcnfg.log" -e "Logs/L$L-b$beta-s%j-g$startcnfg.log" --time=$runtime --mem=$memory submit.script $filepath/run $datpath/GRT $startcnfg $endcnfg $inc $Nc $beta $L $T $xiR $multigen $nTherm $gFixTol $measure $Rmax $Tmax

else

for i in $(seq $startcnfg $(($meas_per_job*$inc)) $endcnfg)
do
endindx=$(($i+($meas_per_job-1)*$inc))
if [ "$endindx" -gt "$endcnfg" ]; then
endindx=$endcnfg
fi

make && sbatch -o "Logs/L$L-b$beta-s%j-m$i.log" -e "Logs/L$L-b$beta-s%j-m$i.log" --time=$runtime --mem=$memory submit.script $filepath/run $datpath/GRT"$i" $i $endindx $inc $Nc $beta $L $T $xiR $multigen $nTherm $gFixTol $measure $Rmax $Tmax

if [ "$endindx" -eq "$endcnfg" ]; then
break
fi

done

fi
