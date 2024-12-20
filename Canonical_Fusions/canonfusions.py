#!/usr/bin/env python
import argparse
import os
import re
import itertools

#./canonfusions.py -input /projects/bgmp/shared/groups/2024/aml/shared/llcombo10_output_files -output /projects/bgmp/shared/groups/2024/aml/lwil/canon_fusions -fusionsfile /projects/bgmp/shared/groups/2024/aml/lwil/BGMP_Student_Project_2024/Canonical_Fusions/canonical_paml_fusions.tsv
  
#Functions
def overlaps(span1, span2):
    return span1[0] <= span2[1] and span2[0] <= span1[1]

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-input", help="Location of folder containing input files", required=True, type=str)
    parser.add_argument("-output", help="Location of folder where the output file will go", required=True, type=str)
    parser.add_argument("-fusionsfile", help="TSV file of canonical fusions", required=True, type=str)
    return parser.parse_args() 

args = get_args()
input = args.input
output = args.output
fusionsfile = args.fusionsfile

#get all file names from folder and save as list
dir_list = os.listdir(input)
patient_list = []

#pull patient IDs & caller names from file names: 
for file in dir_list:
    x = re.split("\\.", file)
    if x[0] not in patient_list: 
        patient_list.append(x[0]) 

found_fusions_dict = {}
potential_fusions_list = []
with open(f"{fusionsfile}", "r") as fusionsfile:
    while True:
        line = fusionsfile.readline()
        line_as_list = line.strip().split("\t")
        
        if line_as_list == ['']:
            break

        key = (line_as_list[0], line_as_list[1], line_as_list[2], line_as_list[3], line_as_list[4])
        potential_fusions_list.append(key)
        #     (gene, mate, chrom): start, end)
        #ex: (CBFA2T3, GLIS2): (chr16, 88941267, 89043615)


#open all files by patient ID 
for patient in patient_list:
    with open(f"{input}/{patient}.vcf", "r") as filein:
        while True:
            line = filein.readline()
            line_as_list = line.strip().split("\t")

            if line_as_list == ['']:
                break

            if line_as_list[2] == "pbsv.BND" or line_as_list[2] == "sniffles.BND":
                bnd_mate_pos = re.search(r"(?<=:)\d+", line_as_list[4])
                bnd_mate_chrom = re.search(r"\w+(?=:)", line_as_list[4])
                temp = (line_as_list[0], line_as_list[1], bnd_mate_chrom.group(), bnd_mate_pos.group())
                #print(patient, value)

                #check first gene
                for fusion in potential_fusions_list:
                    if fusion[2] == temp[0] and temp[1] >= fusion[3] and temp[1] <= fusion[4]:
                        if fusion[1] == "X":
                            print("gene 1 matches", patient, temp, fusion[0], fusion[1])
                        else: 
                            print("gene 1 matches, checking gene 2")
                            matched_gene = fusion[1]
                            #check matched gene
                            for fusion2 in potential_fusions_list:
                                if fusion2[0] == matched_gene and temp[2] == fusion2[2] and temp[3] >= fusion2[3] and temp[3] <= fusion2[4]:
                                    print("gene 2 matches", patient, temp, fusion[0], fusion2[0])
                    
                    elif fusion[2] == temp[2] and temp[3] >= fusion[3] and temp[3] <= fusion[4]:
                        if fusion[1] == "X":
                            print("gene 1 matches", patient, temp, fusion[0], fusion[1])
                        else: 
                            print("gene 1 matches, checking gene 2")
                            matched_gene = fusion[1]
                            #check matched gene
                            for fusion2 in potential_fusions_list:
                                if fusion2[0] == matched_gene and temp[0] == fusion2[2] and temp[1] >= fusion2[3] and temp[1] <= fusion2[4]:
                                    print("gene 2 matches", patient, temp, fusion[0], fusion2[0])




#with open(f"{output}/{patient}_CONSENSUS.vcf", "w") as fileout:
