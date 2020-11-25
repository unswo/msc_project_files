#!/bin/bash
# Open-cravat
# Visualise variants

#Install modules 
oc module install-base
oc module install cosmic
oc module install clinvar
# creates sqlite database from VCF files and prioritises them in accordance with installed modules
oc run $(cat vcfs.txt) -l hg38

oc gui ${NAME_OF_DB}sqlite
