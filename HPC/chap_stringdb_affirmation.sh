#! /bin/bash

#$ -q shai.q
#$ -cwd
#$ -N IL
#$ -l h_vmem=2G

LD_LIBRARY_PATH=/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/lib/:$LD_LIBRARY_PATH

Rscript chap_stringdb_affirmation.r $1
# $1 := the index in the list of chaps to run. (values: 1-15)

