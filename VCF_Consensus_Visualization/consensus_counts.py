#!/usr/bin/env python
import argparse
import os
import re

# ./consensus_counts.py -f /projects/bgmp/shared/groups/2024/aml/shared/llcombo8_output_files -o {your output}

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", help="Location of folder containing input files", required=True, type=str)
    parser.add_argument("-o", help="output counts file name", required=True, type=str)
    return parser.parse_args() 

args = get_args()
f = args.f
o = args.o

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
for cf in dir_list: 
    with open(f"{f}/{cf}", "r") as file: 
        patient = cf.split("_")
        patient = patient[0]
        print(patient)
        while True: 
            line = file.readline().strip().split("\t")

            if line == ['']:
                break 

            if line[0].startswith("##") != True:
                chromosome = line[0]
                #this caller line will need to be changed when we reformat the caller column to be caller.type
                caller = line[2]
                #count the number of svs for each chromosome & add to dictionary
                for sv in svlist: 
                    if line[7].startswith(f"SVTYPE={sv}") or line[7].startswith(f"PRECISE;SVTYPE={sv}") or line[7].startswith(f"IMPRECISE;SVTYPE={sv}"):
                        if (patient, caller, chromosome, sv) in counter_dict:
                            counter_dict[(patient, caller, chromosome, sv)] += 1
                        else: 
                            counter_dict[(patient, caller, chromosome, sv)] = 1 

# write to output file for ALL patients 
with open(o, "w") as counts: 
    counts.write(f"PatientID\tCaller\tChromosome\tSV\tCount\n")
    for key in counter_dict: 
        counts.write(f'{key[0]}\t{key[1]}\t{key[2]}\t{key[3]}\t{counter_dict[key]}\n')