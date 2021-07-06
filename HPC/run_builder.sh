#! /bin/bash

#$ -q shai.q
#$ -cwd
#$ -N IL
#$ -l h_vmem=2G


Rscript build_binari_table.R $1
# $1 := the index in the list of cancers to run. (values: 1-12)


