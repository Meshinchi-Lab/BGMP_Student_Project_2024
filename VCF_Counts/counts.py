#!/usr/bin/env python
import argparse
import os
import re

# ./counts.py -f /projects/bgmp/shared/groups/2024/aml/lwil/VCF_Trial

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", help="Location of folder containing input files", required=True, type=str)
    return parser.parse_args() 

args = get_args()
f = args.f

#get all file names from folder and save as list
dir_list = os.listdir(f)

#initalize
patient_list = []
caller_list = []

#pull patient IDs & caller names from file names: 
for file in dir_list:
    x = re.split("\\.", file)
    if x[0] not in patient_list: 
        patient_list.append(x[0]) 
    if x[1] not in caller_list: 
        caller_list.append(x[1])

#initalize
svlist = ("INS", "DEL", "DUP", "BND", "INV")
counter_dict = dict()

#open each fil by patient ID & caller
for patient in patient_list:
    for caller in caller_list: 
        with open(f"{f}/{patient}.{caller}.vcf", "r") as file: 
            while True: 
                line = file.readline().strip().split("\t")

                if line == ['']:
                    break 
                
                if line[0].startswith("##") != True:
                    chromosome = line[0]

                    #count the number of svs for each chromosome & add to dictionary
                    for sv in svlist: 
                        if line[7].startswith(f"SVTYPE={sv}") or line[7].startswith(f"PRECISE;SVTYPE={sv}") or line[7].startswith(f"IMPRECISE;SVTYPE={sv}"):
                            if (caller, chromosome, sv) in counter_dict:
                                counter_dict[(caller, chromosome, sv)] += 1
                            else: 
                                counter_dict[(caller, chromosome, sv)] = 1 

#write to output file for graphing
with open("counts.tsv", "w") as counts: 
    counts.write("Caller\tChromosome\tSV\tCount\n")
    for key in counter_dict: 
        counts.write(f'{key[0]}\t{key[1]}\t{key[2]}\t{counter_dict[key]}\n')