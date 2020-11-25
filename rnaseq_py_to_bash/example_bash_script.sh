#!/bin/bash 
#SBATCH --array=0-6 
#SBATCH --cores=10
#SBATCH --nodes=1
#SBATCH --mem=60g
#SBATCH --partition=defq

declare -a read_array=("SRRexamplefile_KO" "SRCexamplefile_WT" "SR3examplefile_KO" "SRDexamplefile_WT" "sample1" "sample2" "sample4")

trim_galore -j 10 -q 30 --basename ${read_array[$SLURM_ARRAY_TASK_ID]} --paired --fastqc_args "--outdir fastqc_output" ${read_array[$SLURM_ARRAY_TASK_ID]}*.fastq.gz 

STAR --genomeDir star_index --readFilesIn ${read_array[$SLURM_ARRAY_TASK_ID]}_val_1.fq.gz ${read_array[$SLURM_ARRAY_TASK_ID]}_val_2.fq.gz --outReadsUnmapped None --runThreadN 10 --outTmpDir $TMPDIR --outFileNamePrefix star_output/${read_array[$SLURM_ARRAY_TASK_ID]}_ --twopassMode Basic --readFilesCommand "gunzip -c" --outSAMstrandField intronMotif --outSAMunmapped Within --chimSegmentMin 12 --chimJunctionOverhangMin 8 --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --outSAMattrRGline ID:GRPundef --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --alignInsertionFlush Right --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 

samtools sort -@ 10 star_output/${read_array[$SLURM_ARRAY_TASK_ID]}_Aligned.out.sam -o ${read_array[$SLURM_ARRAY_TASK_ID]}.sorted.bam  
 mv ${read_array[$SLURM_ARRAY_TASK_ID]}.sorted.bam star_output/sorted_bams/

STAR-Fusion --genome_lib_dir ctat_fusion_lib -J star_output/${read_array[$SLURM_ARRAY_TASK_ID]}_Chimeric.out.junction --output_dir fusion_output/${read_array[$SLURM_ARRAY_TASK_ID]} 

FusionInspector --fusions fusion_output/${read_array[$SLURM_ARRAY_TASK_ID]}/star-fusion.fusion_predictions.abridged.tsv  --out_prefix ${read_array[$SLURM_ARRAY_TASK_ID]} --min_junction_reads 1  --min_novel_junction_support 3 --min_spanning_frags_only 5 --vis --max_promiscuity 10 --output_dir fusion_output/${read_array[$SLURM_ARRAY_TASK_ID]}/FusionInspector --genome_lib_dir ctat_fusion_lib --CPU 10 --include_Trinity  --annotate --left_fq reads_dir/${read_array[$SLURM_ARRAY_TASK_ID]}_val_1.fq.gz --right_fq reads_dir/${read_array[$SLURM_ARRAY_TASK_ID]}_val_2.fq.gz