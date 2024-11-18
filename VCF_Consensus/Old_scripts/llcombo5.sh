#!/bin/bash

#SBATCH --partition=bgmp
#SBATCH --account=bgmp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lenara@uoregon.edu
#SBATCH --output=llcombo_%j.out

/usr/bin/time -v ./llcombo5.py -i /projects/bgmp/shared/groups/2024/aml/lenara/trial_files -o /projects/bgmp/shared/groups/2024/aml/lenara/BGMP_Student_Project_2024/VCF_Consensus

