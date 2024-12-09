#!/bin/bash
# $1 cnfg filepath
# $2 data filepath
# $3 cnfg start num
# $4 cnfg end num
# $5 cnfg increment
# $6 N colors
# $7 beta
# $8 L
# $9 T
# $10 xi_R
# $11 (bool) generating initial cnfgs
# $12 number of thermalization sweeps
# $14 dF for coulomb gauge
# $15 (bool) isMeasuring
# $16 R_max
# $17 T_max
# $18 time limit
# $19 memory req
#makes the folders for the configs and data

cnfgpath=$1/run
datpath=$2/run

sbatch --time=${19} --mem=${18} submit.script $cnfgpath $datpath $3 $4 $5 $6 $7 $ 8 $9 ${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17}


