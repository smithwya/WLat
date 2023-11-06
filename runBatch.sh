#!/bin/bash
# $1 number of jobs to submit
# $2 beta
# $3 N (lattice spatial extent)
# $4 T (lattice temporal extent)
# $5 xiR
# $6 1=hot, 0=cold (continue thermalizing current file or start over)
# $7 number of sweeps to perform
# $8 tolerance of gauge-fixing (dF<= tol)
# $9 suffix (name to add to end of datafiles)
# $10 number of measurements
# $11 sweeps between each measurement
# $12 runtime
# $13 filepath
# $14 Nc

#makes the folders for the configs and data
cd ${13}
mkdir -p ./Configs/"$3"^3x"$4"/$2-$5
mkdir -p ./Data/"$3"^3x"$4"/$2-$5
configpath="${13}"/Configs/"$3"^3x"$4"/$2-$5
datapath="${13}"/Data/"$3"^3x"$4"/$2-$5

cd -

#submits jobs
for ((i=1; i<= $1; i++))
do
#submit to carbonate:
sbatch --time=${12} submit.script "$i" $2 $3 $4 $5 $6 $7 $8 $configpath $datapath $9 ${10} ${11} ${14}

#for testing
#./bin/WLat "$i" $2 $3 $4 $5 $6 $7 $8 $configpath $datapath $9 ${10} ${11} ${14}
done
