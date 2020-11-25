#!/bin/bash
#
# variant calling
#SBATCH -c 10
#SBATCH --mem 60g
#SBATCH --mail-type=ALL
#

source ~/.bashrc

conda activate CTAT_mut

ctat_mutations --left ${1}_val_1.fq.gz \
--right ${1}_val_2.fq.gz \
--out_dir $PROJECT/var_output \
--genome_lib_dir $PROJECT/CTAT_lib/GRCh38_gencode_v33_CTAT_lib_Apr062020.source/ctat_genome_lib_build_dir \
--threads 10 \
--star_memory 40g

conda deactivate

