#!/usr/bin/env python
# Adrienne Unsworth, 2020
from pipeline_functions import *
from tools import *

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
