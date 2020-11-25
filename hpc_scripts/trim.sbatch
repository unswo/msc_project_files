#!/bin/bash
#
# trim with trimgalore
#SBATCH -c 4
#SBATCH --mem 5g
#SBATCH --mail-type=END

module load Java
module load cutadapt

trim_galore -j 4 -q 20 --basename $(basename ${1}) --paired --fastqc_args "--outdir $PROJECT/data/trim_galore" $(basename ${1})_*_1.fastq.gz $(basename ${1})_*_2.fastq.gz

module purge 
