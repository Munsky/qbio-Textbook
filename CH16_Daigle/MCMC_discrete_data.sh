#!/bin/bash

####
#### Estimates parameters from incomplete data using MCMC method
####

if [ $# -lt 6 ]; then
	echo "Usage: `basename $0` SBMLFILE(character) DATAFILE(character) NUMITER(numeric) NUMBURN(numeric) NUMTHIN(numeric) TABLEFILE(character) [SEED](numeric)"
	exit 65
fi
sbmlfile=$1
datafile=$2
numiter=$3
numburn=$4
numthin=$5
tablefile=$6
if [ $# -ge 7 ]; then
	seed=$7
else
	seed=101
fi

### Store executable paths
SRC=`dirname $0`
UTILS="${SRC}/utils"

### Run MCMC; burn-in and thin output
numobs=`tail -n +2 ${datafile} | wc -l`
${SRC}/stochInf/src/stochInf -m ${sbmlfile} -d ${datafile} -n ${numobs} -i ${numiter} -t .05 -s ${seed} -v 0 | ${UTILS}/mcmc.py -b ${numburn} -t ${numthin} > ${tablefile}

### Rescale parameters based on data time points
numparams=`grep -c "parameter" ${sbmlfile}`
finaliter=`tail -n1 ${datafile} | cut -f1 -d' '`
finaltime=`tail -n1 ${datafile} | cut -f2 -d' '`
multiplier=`echo "${finaliter}/${finaltime}" | bc -l`
${SRC}/scale_MCMC_parameters.R ${tablefile} ${numparams} ${multiplier}

### Remove glpk files
rm glpk.log glpk.txt

