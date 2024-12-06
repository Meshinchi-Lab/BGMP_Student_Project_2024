#!/bin/bash

#SBATCH --partition=bgmp
#SBATCH --account=bgmp
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lenara@uoregon.edu

/usr/bin/time -v ./consensus_counts.py -f /projects/bgmp/shared/groups/2024/aml/shared/llcombo10_output_files -o llcombo10_consensus_counts.tsv