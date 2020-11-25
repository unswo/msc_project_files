#!/bin/bash
#
# pipeline
#SBATCH -c 10
#SBATCH --mem 60g
#SBATCH --mail-type=ALL
#SBATCH--array=0-62
# Assumes all of these tools are availble on the command line or using the conda environments provided

# work around for samtools dependencies
# change into the lib directory for the conda environment
# cd ~/.conda/envs/
# ln -s libcrypto.so.1.1 libcrypto.so.1.0.0

# array=0-61:2

declare -a normals=("SRR6059542" "SRR6059543" "SRR6059546" "SRR6059547" "SRR6059550" "SRR6059551" "SRR6059552" "SRR6059553" "SRR6059558" "SRR6059559" "SRR6059562" "SRR6059563" "SRR6059566" "SRR6059567" "SRR6059568" "SRR6059569" "SRR6059574" "SRR6059575" "SRR6059578" "SRR6059579" "SRR6059582" "SRR6059583" "SRR6059584" "SRR6059585" "SRR6059590" "SRR6059591" "SRR6059592" "SRR6059593" "SRR6059596" "SRR6059597" "SRR6059602" "SRR6059603" "SRR9623678" "SRR9623679" "SRR9623682" "SRR9623683" "SRR9623686" "SRR9623687" "SRR9623690" "SRR9623691" "SRR9623694" "SRR9623695" "SRR9623698" "SRR9623699" "SRR9623702" "SRR9623703" "SRR9623706" "SRR9623707" "SRR9623710" "SRR9623711" "SRR9623714" "SRR9623715" "SRR9623718" "SRR9623719" "SRR9623722" "SRR9623723" "SRR9623726" "SRR9623727" "SRR9623730" "SRR9623731" "SRR9623734" "SRR9623735")

declare -a tumours=("SRR6059540" "SRR6059541" "SRR6059544" "SRR6059545" "SRR6059548" "SRR6059549" "SRR6059554" "SRR6059555" "SRR6059556" "SRR6059557" "SRR6059560" "SRR6059561" "SRR6059564" "SRR6059565" "SRR6059570" "SRR6059571" "SRR6059572" "SRR6059573" "SRR6059576" "SRR6059577" "SRR6059580" "SRR6059581" "SRR6059586" "SRR6059587" "SRR6059588" "SRR6059589" "SRR6059594" "SRR6059595" "SRR6059598" "SRR6059599" "SRR6059600" "SRR6059601" "SRR9623676" "SRR9623677" "SRR9623680" "SRR9623681" "SRR9623684" "SRR9623685" "SRR9623688" "SRR9623689" "SRR9623692" "SRR9623693" "SRR9623696" "SRR9623697" "SRR9623700" "SRR9623701" "SRR9623704" "SRR9623705" "SRR9623708" "SRR9623709" "SRR9623712" "SRR9623713" "SRR9623716" "SRR9623717" "SRR9623720" "SRR9623721" "SRR9623724" "SRR9623725" "SRR9623728" "SRR9623729" "SRR9623732" "SRR9623733")

source ~/.bashrc

conda activate variant_calling

samtools sort -@ 10 ${tumours[$SLURM_ARRAY_TASK_ID]}_Aligned.out.sam -o ${tumours[$SLURM_ARRAY_TASK_ID]}.sorted.bam

picard MarkDuplicates \
I=${tumours[$SLURM_ARRAY_TASK_ID]}.sorted.bam \
O=${tumours[$SLURM_ARRAY_TASK_ID]}_marked_dup.bam \
M=${tumours[$SLURM_ARRAY_TASK_ID]}_marked_duo_metrics.txt

gatk --java-options "-Xmx50G" SplitNCigarReads -I ${tumours[$SLURM_ARRAY_TASK_ID]}_marked_dup.bam -O ${tumours[$SLURM_ARRAY_TASK_ID]}.split_sort.bam -R $FASTA

gatk --java-options "-Xmx50G" BaseRecalibrator -I ${tumours[$SLURM_ARRAY_TASK_ID]}.split_sort.bam -R $PROJECT/CTAT_lib/GRCh38_gencode_v33_CTAT_lib_Apr062020.source/ctat_genome_lib_build_dir/ref_genome.fa --known-sites $COSMIC_VCF --known-sites $DBSNP_VCF -O ${tumours[$SLURM_ARRAY_TASK_ID]}_recal.table

gatk --java-options "-Xmx50G" ApplyBQSR --bqsr-recal-file ${tumours[$SLURM_ARRAY_TASK_ID]}_recal.table --input ${tumours[$SLURM_ARRAY_TASK_ID]}.split_sort.bam --output ${tumours[$SLURM_ARRAY_TASK_ID]}.recal.bam

gatk --java-options "-Xmx50G" Mutect2 -R $FASTA \
-I ${tumours[$SLURM_ARRAY_TASK_ID]}.recal.bam \
-I ${tumours[$(($SLURM_ARRAY_TASK_ID + 1))]}.recal.bam \
-I ${normals[$SLURM_ARRAY_TASK_ID]}.recal.bam \
-I ${normals[$(($SLURM_ARRAY_TASK_ID + 1))]}.recal.bam \
-normal ${normals[$SLURM_ARRAY_TASK_ID]} \
-normal ${normals[$(($SLURM_ARRAY_TASK_ID + 1))]} \
--germline-resource $GNOMAD_VCF \
--panel-of-normals pon.vcf.gz \
-O ${tumours[$SLURM_ARRAY_TASK_ID]}_${tumours[$(($SLURM_ARRAY_TASK_ID + 1))]}.vcf.gz \
--tmp-dir $TMPDIR

gatk --java-options "-Xmx50G" IndexFeatureFile -I ${tumours[$SLURM_ARRAY_TASK_ID]}_${tumours[$(($SLURM_ARRAY_TASK_ID + 1))]}.vcf.gz

gatk --java-options "-Xmx50G" FilterMutectCalls -R $FASTA \
-V ${tumours[$SLURM_ARRAY_TASK_ID]}_${tumours[$(($SLURM_ARRAY_TASK_ID + 1))]}.vcf.gz \
-O ${tumours[$SLURM_ARRAY_TASK_ID]}_${tumours[$(($SLURM_ARRAY_TASK_ID + 1))]}_filtered.vcf.gz \
--tmp-dir $TMPDIR

bcftools index -f ${tumours[$SLURM_ARRAY_TASK_ID]}_${tumours[$(($SLURM_ARRAY_TASK_ID + 1))]}_filtered.vcf.gz
bcftools filter ${tumours[$SLURM_ARRAY_TASK_ID]}_${tumours[$(($SLURM_ARRAY_TASK_ID + 1))]}_filtered.vcf.gz -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chrX,chrY > ${tumours[$SLURM_ARRAY_TASK_ID]}_${tumours[$(($SLURM_ARRAY_TASK_ID + 1))]}.processed.vcf
