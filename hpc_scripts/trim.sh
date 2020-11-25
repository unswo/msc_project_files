#!/bin/bash
#launch trimgalore extra wastefully

for x in $PROJECT/data/PRJNA411786/SRR*;
do
cd $x
z=$(basename $x)
sbatch ~/scripts/trim.sbatch $z
cd ..
done
 
