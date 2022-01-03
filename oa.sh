#!/bin/bash

if [ $# -lt 1 ];then
	echo "Usage: sh $0 <sample_ID> <reads1> <read2>"
	exit
fi

oa --log $1\.index.log index --length 125 $1 $2 $3
oa --log $1\.buildgraph.log buildgraph --minread 5 --assmax 200000 --smallbranches 15 --seeds protChloroArabidopsis $1 $1\.chlo
oa --log $1\.unfold.log unfold $1 $1\.chlo 1>$1\.chlo.fa
