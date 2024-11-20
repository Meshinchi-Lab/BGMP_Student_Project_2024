#!/bin/bash

#SBATCH --partition=bgmp
#SBATCH --account=bgmp
#SBATCH --mem=160G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lenara@uoregon.edu

/usr/bin/time -v ./consensus_counts.py -f /projects/bgmp/shared/groups/2024/aml/shared/llcombo8_output_files -o llcombo8_consensus_counts.tsv