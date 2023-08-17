#!/bin/bash
# $1 number of jobs to submit
# $2 beta
# $3 N (lattice spatial extent)
# $4 T (lattice temporal extent)
# $5 xiR
# $6 1=hot, 0=cold (continue thermalizing current file or start over)
# $7 number of sweeps to perform
# $8 tolerance of gauge-fixing (dF<= tol)

#makes the folders for the configs and data
cd /N/slate/smithwya/"SU(2)"
mkdir -p ./Configs/"$3"^3x"$4"/$2-$5
mkdir -p ./Data/"$3"^3x"$4"/$2-$5
configpath=/N/slate/smithwya/"SU(2)"/Configs/"$3"^3x"$4"/$2-$5
datapath=/N/slate/smithwya/"SU(2)"/Data/"$3"^3x"$4"/$2-$5

cd -

#compiles code
make

#submits jobs
for ((i=1; i<= $1; i++))
do
#submit to carbonate:
sbatch submit.script "$i" $2 $3 $4 $5 $6 $7 $8 $configpath $datapath

#for testing: 
#./bin/WLat "$i" $2 $3 $4 $5 $6 $7 $8 $configpath $datapath
done