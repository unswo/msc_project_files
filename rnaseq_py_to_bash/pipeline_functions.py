#!/usr/bin/env python
# Adrienne Unsworth, 2020
import os
import configparser
import subprocess
import re


def load_config(configfile='config.cfg'):
    '''Load config file'''
    config = configparser.ConfigParser()
    config.read(configfile)
    return config


def load_samples(configfile='config.cfg'):
    '''Load samples'''
    config = load_config(configfile)
    samples = dict(config.items('SAMPLES'))
    sample_name_file = samples['sample_names']
    sample_names = open(sample_name_file)
    sample_names = sample_names.readlines()
    sample_names = [x.strip('\n') for x in sample_names]
    return sample_names


def list_to_bash_array(array_name='read_array', configfile='config.cfg'):
    '''Convert text file of samples separated by new lines to a bash array'''
    sample_names = load_samples()
    num_samples = len(sample_names)
    print('{0} samples in file.'.format(num_samples))
    sample_names = '"' + '" "'.join(sample_names) + '"'
    cmd = 'declare -a {0}=(' + sample_names + ')'
    return cmd.format(array_name), num_samples


def getdir(configfile='config.cfg'):
    '''Get relevant directories relating to tools, reads and output.'''
    config = load_config(configfile)
    directories = dict(config.items('DIRECTORIES'))
    for i in directories:
        directories[i] = directories[i].split('#', 1)[0].strip()
        globals().update(directories)
    print('Checking resource paths exist...')
    if not os.path.exists(reads_dir):
        print('Read directory not found. Please correct in the config file.')
        exit(1)
    if not os.path.exists(ctat_fusion_lib):
        print(
            'CTAT fusion library not found. Please check input directory or generate library using the available script.')
        exit(1)
    if not os.path.exists(star_index):
        print('Star index not found')
        exit(1)
    print('Ok!')
    print('Creating output directories if they dont exist...')
    os.system('mkdir -p {0}/sorted_bams {1} {2} {3}'.format(star_output, fastqc_output, fusion_output, var_output))
    return directories


def sbatch_settings(configfile='config.cfg'):
    '''Load SLURM settings for sbatch command from config file.'''
    config = load_config(configfile)
    sample_list, num_sample = list_to_bash_array()
    s_settings = dict(config.items('SBATCH'))
    array_size = num_sample - 1
    cmd ='#!/bin/bash \n'
    if s_settings['array'] == 'default':
        cmd = cmd + '#SBATCH --array=0-{0} \n' .format(array_size)
        del s_settings['array']
    for i in s_settings.keys():
        s_settings[i] = s_settings[i].split('#', 1)[0].strip()
        cmd = cmd + '#SBATCH --' + i + '=' + s_settings[i] + '\n'

    return s_settings, cmd


def pipesettings(configfile='config.cfg'):
    '''Load pipeline settings such as what tools to use and whether we have paired end data.'''
    config = load_config(configfile)
    p_settings = config['PIPELINE']
    print('Using the following settings:')
    for x in p_settings.keys():
        if p_settings.getboolean(x):
            print(x)
    paired_end = p_settings.getboolean('paired_end')
    return p_settings, paired_end



