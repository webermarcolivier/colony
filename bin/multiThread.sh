#!/bin/bash

# Script file for running multiple simulations in batch on a local computer.
# This script lists all the input files in the current directory,
# finds the executable in the current directory, and runs a simulation
# with the same executable for each input file.
# The number of running jobs are controlled with the GNU parallel package.
# (see below).
# Author: Marc Weber
# Last update: 2015.08.21

NB_CPU_MAX=4
MEMORY_MAX=500M

# Count the number of input files
INPUTFILES=`find . -maxdepth 1 -name "input???.dat" -type f | sort`
COUNT=`find . -maxdepth 1 -name "input???.dat" -type f | wc -l`
echo
echo "Number of input files found: " $COUNT

# Extract the id number of all the input files and store it in an array
IFILE=1
for f in $INPUTFILES; do
	IDLIST[$IFILE]=`expr match "$f" '.*\([0-9][0-9][0-9]\)'`
	echo ${IDLIST[$IFILE]}
	let IFILE=$IFILE+1
done
IDLISTSTRING=$( IFS=$'\n'; echo "${IDLIST[*]}" )

# Find the executable
COUNT=`find . -maxdepth 1 -name "*.x" -type f | wc -l`
if [ $COUNT -ne 1 ]; then
	echo "ERROR: There is less than or more than 1 executable file *.x in the current directory."
	exit 0 
fi
EXECUTABLE=`ls *.x`
echo "Executable file: " $EXECUTABLE

# Launch each thread with the same executable and a different input file
# Use the programm GNU parallel to limit the number of cores used for the
# simulation batch. 
# More information: http://www.gnu.org/software/parallel
# option --xapply do a one-to-one mapping between the two input sources
# option -j sets the maximum number of cores
# option --memfree size minimum memory free when starting another job
parallel -vv --xapply -j $NB_CPU_MAX --memfree $MEMORY_MAX ./run.sh $EXECUTABLE {1} {2} ::: $INPUTFILES ::: $IDLISTSTRING


# Otherwise: launch all processes at the same time
#ITHREAD=1
#for f in $INPUTFILES; do
#	INPUTFILENAME=$(basename $f)
#	echo "Launching simulation, thread #""$ITHREAD"" with input file: ""$INPUTFILENAME"
#	nice nohup ./"$EXECUTABLE" "$INPUTFILENAME" > out"$ITHREAD".txt & 
#	let ITHREAD=ITHREAD+1
#done

wait

