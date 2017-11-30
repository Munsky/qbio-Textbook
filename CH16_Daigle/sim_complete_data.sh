#!/bin/bash

####
#### Simulates complete data trajectory using gillespie executable
####

if [ $# -lt 3 ]; then
	echo "Usage: `basename $0` SBMLFILE(character) TIME(numeric) OUTFILE(character) [SEED](numeric)"
	exit 65
fi
sbmlfile=$1
time=$2
outfile=$3
if [ $# -ge 4 ]; then
	seed=$4
else
	seed=101
fi

### Simulate complete data trajectory
SRC=`dirname $0`
${SRC}/gillespie/src/gillespie -m ${sbmlfile} -t ${time} -s ${seed} | sed 's/ $//' > ${outfile}

### Append final system state if needed
tf=`tail -n1 ${outfile} | cut -f2 -d' '`
finished=`echo "${tf} == ${time}" | bc`
if [ ${finished} -eq 0 ]; then
	iter=`tail -n +2 ${outfile} | wc -l`
	final=`mktemp`
	tail -n1 ${outfile} > ${final}
	sed -i -e "s/^[0-9]* [0-9.]* \(.*\)$/${iter} ${time} \1/" ${final}
	echo "`cat ${outfile} ${final}`" > ${outfile}
fi

