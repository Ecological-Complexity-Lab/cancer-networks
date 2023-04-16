#! /bin/bash

#$ -q shai.q
#$ -cwd
#$ -N IL

LD_LIBRARY_PATH=/gpfs0/shai/projects/R4/R-4.2.0/lib64/R/lib/:$LD_LIBRARY_PATH

Rscript chap_stringdb_affirmation.r $1
# $1 := the index in the list of chaps to run. (values: 1-15)

