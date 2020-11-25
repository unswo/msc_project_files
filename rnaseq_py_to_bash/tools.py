#!/usr/bin/env python
# Adrienne Unsworth, 2020
from pipeline_functions import *


def trimgalore(cores, paired_end, fastqc_output, read='${read_array[$SLURM_ARRAY_TASK_ID]}'):
    cmd = 'trim_galore -j {0} -q 30 --basename {1} {2} ' \
          '--fastqc_args "--outdir {3}" {1}*.fastq.gz '
    if paired_end:
        return cmd.format(cores, read, '--paired', fastqc_output)
    else:
        return cmd.format(cores, read, '', fastqc_output)


def star(paired_end, star_index, star_output, cores, fusion, two_pass, rsem_use, configfile='config.cfg',
         read='${read_array[$SLURM_ARRAY_TASK_ID]}'):
    config = load_config(configfile)
    tmpdir = config['ENVIRONMENT']['tmp_dir']
    sample = '{0}_val_1.fq.gz {0}_val_2.fq.gz'
    twopass = ''
    rsem = ''
    if rsem_use:
        rsem = '--quantMode TranscriptomeSAM'
    if fusion:
        # uses star variables as defined at the star fusion wiki on github
        if not paired_end:
            sample = '{0}.trimmed.fq.gz'
        sample = sample.format(read)
        cmd = 'STAR --genomeDir {0} --readFilesIn {1} ' \
              '--outReadsUnmapped None --runThreadN {2} --outTmpDir {3} --outFileNamePrefix {4}/{5}_ ' \
              '--twopassMode Basic --readFilesCommand "gunzip -c" --outSAMstrandField intronMotif ' \
              '--outSAMunmapped Within --chimSegmentMin 12 --chimJunctionOverhangMin 8 --chimOutJunctionFormat 1 ' \
              '--alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 ' \
              '--outSAMattrRGline ID:GRPundef --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 ' \
              '--chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --alignInsertionFlush Right --alignSplicedMateMapLminOverLmate 0 ' \
              '--alignSplicedMateMapLmin 30 ' \
              '{6}'.format(star_index, sample, cores, tmpdir, star_output, read, rsem)
        return cmd
    if two_pass:
        twopass = '--twopassMode Basic'
    if not paired_end:
        sample = '{0}.trimmed.fq.gz'
    sample = sample.format(read)
    cmd = 'STAR --genomeDir {0} --readFilesIn {1} ' \
          '--outReadsUnmapped None --runThreadN {2} --outTmpDir {3} --outFileNamePrefix {4}/{6}_ ' \
          '{5} --readFilesCommand "gunzip -c" {7}'.format(star_index, sample, cores, tmpdir, star_output, twopass, read, rsem)
    return cmd


def star_fusion(ctat_fusion_lib, reads_dir, star_output, fusion_output, fusion_inspector, cores,
                read='${read_array[$SLURM_ARRAY_TASK_ID]}'):
    cmd = 'STAR-Fusion --genome_lib_dir {0} -J {1}/{2}_Chimeric.out.junction --output_dir {3}/{2} ' \
        .format(ctat_fusion_lib, star_output, read, fusion_output, )
    if fusion_inspector:
        fus_cmd = 'FusionInspector --fusions {0}/{1}/star-fusion.fusion_predictions.abridged.tsv  ' \
                  '--out_prefix {1} --min_junction_reads 1  --min_novel_junction_support 3 --min_spanning_frags_only 5 ' \
                  '--vis --max_promiscuity 10 --output_dir {0}/{1}/FusionInspector ' \
                  '--genome_lib_dir {2} --CPU {3} --include_Trinity  --annotate ' \
                  '--left_fq {4}/{1}_val_1.fq.gz --right_fq {4}/{1}_val_2.fq.gz' \
            .format(fusion_output, read, ctat_fusion_lib, cores, reads_dir)
        cmd = cmd + '\n' * 2 + fus_cmd
    return cmd


def samtools_sort(star_output, cores, read='${read_array[$SLURM_ARRAY_TASK_ID]}'):
    cmd = 'samtools sort -@ {2} {1}/{0}_Aligned.out.sam -o {0}.sorted.bam  \n ' \
          'mv {0}.sorted.bam {1}/sorted_bams/' \
        .format(read, star_output, cores)
    return cmd

# pull settings from config file
directories = getdir()
globals().update(directories)
p_settings, paired_end = pipesettings()
samples_list = load_samples()
environ_settings, sbatch_args = sbatch_settings()

# declare empty commands
array_cmd, num_samples = list_to_bash_array()
trim_cmd = ''
star_cmd = ''
fusion_cmd = ''

# build script
if p_settings.getboolean('qc_trimgalore'):
    trim_cmd = trimgalore(paired_end = paired_end, cores = environ_settings['cores'], fastqc_output = fastqc_output)
if p_settings.getboolean('star'):
    star_cmd = star(paired_end, star_index, star_output, environ_settings['cores'], rsem_use = p_settings.getboolean('rsem'), fusion = p_settings.getboolean('star_fusion'), two_pass = p_settings.getboolean('2_pass'))
    star_cmd = star_cmd + '\n'*2 + samtools_sort(star_output, environ_settings['cores'])
    if p_settings.getboolean('star_fusion'):
        fusion_cmd = star_fusion(ctat_fusion_lib, reads_dir, star_output, fusion_output, p_settings.getboolean('fusion_inspector'), environ_settings['cores'])


# write script
bash_script = sbatch_args + '\n' + array_cmd + '\n'*2 +  trim_cmd + '\n'*2 + star_cmd + '\n'*2 + fusion_cmd
with open('example_bash_script.sh', 'w') as file:
    file.write(bash_script)
