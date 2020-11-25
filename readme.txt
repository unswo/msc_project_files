# Source code for 'A meta-analysis of publicly available RNA-Seq data in Prostate Cancer'
# Adrienne Unsworth 190033848, 2020

# Data
Data (raw fastqs) is available at the SRA under BioProject PRJNA411786 and PRJNA552058. 
Details of exact metadata used can be seen under proj_design.csv. fusions.csv contains fusions identified from the worklfow following in silico validation by FusionInspector and variants.csv in the variant folder contains output from the variant calling section of the pipeline. DE_genes.csv and hallmark_GSEA.csv in the DEGs folder contain DE genes identified from the analysis and results from gene set enrichment analysis respectively.
Unfortunately no toy datasets are included to demonstrate how the pipeline works as each transcript file and the results R objects are >5MB. Scripts are purely explanatory and build upon the methods of the dissertation.

# Scripts
All tools are assumed to be available on the command line.
Align_Fusion_AltSplice.sh includes all of the commands used for the project, with the exception of variant calling which is multiple scripts in the 'variants' folder. Does not include commands to generate indexes for STAR and CTAT genome library. More info on each staep is available in the script.

normals_for_PON.sbatch generates vcfs for normal samples and create_pon_db.sbatch compiles a database from all of these VCFs to create a panel of normals. tumour_mutext.sbatch can then be run in order to identify somatic variants. opencravat.sh then collates the VCFs produced into an sqlite database.

#Environments
envs contains 2 conda environments. rnaseqpipe is for trimming, alignment, transcript quantification and fusion detection. variant_calling is for variant calling using the scripts in the variants folder.
No environment is provided for use with MAJIQ as it requires an academic license - see the MAJIQ website for more details. (Note - needs installing into an empty conda environment via pip)
R_session_info.txt contains version and repository information for R packages used for DEGs_workflow.R in the DEGs folder.

A work around is necessary for samtools dependencies: change into the lib directory for the conda environment
Add symlink for needed library:
ln -s libcrypto.so.1.1 libcrypto.so.1.0.0

#rnaseq_py_to_bash
rnaseq_py_to_bash is a set of (incomplete) scripts to generate a batch file for slurm from a config file and a text file containing a list of sample names separated by newlines. This was intended to make the workflow described in my dissertation modular and reproducible for (potential) use with the RNA-seq data from the original planned wet-lab. The scripts were never updated after the switch from using RSEM as a quantification method to Stringtie and therefore does not replicate the workflow exactly nor does it in it's current form have any real error handling. 
Uses configparser to import a cfg file, select tools from those available in tools.py and populates the strings encoding commands for each tools with input from the config file. Not all options have been implemented properly and are commented as such in the example config file.

##Usage
Run main.py after populating config.cfg and samples.txt. Assumes tools are available on the command line a relevant indexes (e.g. STAR index) are available.

