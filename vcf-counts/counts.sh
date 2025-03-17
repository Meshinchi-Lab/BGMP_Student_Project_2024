#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --time=1-0

/usr/bin/time -v ./counts.py -f /projects/bgmp/shared/groups/2024/aml/shared/Our_Variant_Call_Files

#Rscript -e "install.packages('rmarkdown')"
#Rscript -e "rmarkdown::render('counts_visualization.rmd')"

#Rscript -e "install.packages('rmarkdown')"

