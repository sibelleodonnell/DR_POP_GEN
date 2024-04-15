#!/bin/bash
# runs a script supplied by the user on a list of samples
# submitting each job to the job schedule
# this version works on USC's HPC with SLURM
# usage: batch.sh script list jobname

SCRIPT=$1
LIST=$2
NAME=$3

if [ "$#" -ne 3 ]
	then
	echo ""
	echo "Usage: batch.sh script list prefix"
	echo "Where"
	echo "	script:		name of the job submission script you want to run on multiple samples"
	echo "	list:		a list of sample names (directory names)"
	echo "	prefix:		a prefix for your job names"
	echo "For example. 'batch.sh proc.sh samples.list proc'"
	echo ""
	exit
fi

exec < "$LIST"
while read LINE
do
        echo "$LINE"
	cd $LINE
	rm -rf $NAME.$LINE
	sbatch -J $NAME.$LINE --export=NAME=$LINE ../$SCRIPT
	cd ../
done

