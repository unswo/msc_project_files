# rnaseq-pipeline-msc
Use a python script to generate a bash script for SLURM. Assumes all tools are available from the command line. A suitable conda environment can be generated from the .yaml file included, though some amendments may be required.

This script is intended to replicate the main workflow used in this MSc project in a modular fashion; i.e. you can choose whether to include certain steps such as trimming and alignment if already performed. Though this is with the caveat of the STAR function will be expecting the fastq files to be trimmed and have a filename reflective of being trimmed with TrimGalore!. This is easily amended by editing/creating new functions. An example annotated function can be seen below.

# How to use
NOTE: You need to set your GTF, and genome FASTA at minimum to use the default config file. 
The example config file will perform a standard workflow of trimming -> alignment with STAR -> quantification with RSEM using the parameters that you have set in the config file.


