#!/usr/bin/env python
import argparse
import os
import re
import itertools

#./canonfusions.py -input /projects/bgmp/shared/groups/2024/aml/shared/llcombo10_output_files -output /projects/bgmp/shared/groups/2024/aml/lwil/canon_fusions -fusionsfile /projects/bgmp/shared/groups/2024/aml/lwil/BGMP_Student_Project_2024/Canonical_Fusions/full_canonical_fusions_list.tsv

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

#FUNCTIONS
def overlaps(span1, span2):
    return span1[0] <= span2[1] and span2[0] <= span1[1]

#MAKE PATIENT LIST
#get all file names from folder and save as list
dir_list = os.listdir(input)
patient_list = []

#pull patient IDs & caller names from file names: 
for file in dir_list:
    x = re.split("\\.vcf", file)
    if x[0] not in patient_list: 
        patient_list.append(x[0]) 

#MAKE FUSION LIST FROM FUSION FILE
#this will be read in as a dataframe in R 
potential_fusions_list = []
with open(f"{fusionsfile}", "r") as fusionsfile:
    while True:
        line = fusionsfile.readline()
        line_as_list = line.strip().split("\t")
        
        if line_as_list == ['']:
            break

        key = (line_as_list[0], line_as_list[1], line_as_list[2], line_as_list[3], line_as_list[4])
        potential_fusions_list.append(key)
        #     (gene, mate, chrom, start, end)
        #ex: (CBFA2T3, GLIS2, chr16, 88941267, 89043615)

found_fusions_list = []

###############################
# THIS IS WHERE I NEED HELP!! #
###############################

#FIND CANONICAL FUSIONS IN EVERY FILE 
#open all files by patient ID 
for patient in patient_list:    #iterate through the list of patient names
    with open(f"{input}/{patient}.vcf", "r") as filein:
        while True:
            line = filein.readline()
            line_as_list = line.strip().split("\t")

            if line_as_list == ['']:
                break

            if line_as_list[0].startswith("##"): #skip header lines
                continue

            if line_as_list[2] == "pbsv.BND" or line_as_list[2] == "sniffles.BND" or line_as_list[2].startswith("Sniffles2.BND") or line_as_list[2].startswith("pbsv.BND"): #if sv type is BND 
                bnd_mate_pos = re.search(r"(?<=:)\d+", line_as_list[4]) #pull out the mate position from the INFO column
                bnd_mate_chrom = re.search(r"\w+(?=:)", line_as_list[4]) #pull out the mate chromosome from the INFO column
                temp = (line_as_list[0], line_as_list[1], bnd_mate_chrom.group(), bnd_mate_pos.group()) #make a temp variable to hold the chr, position, mate chr, and mate position
                #temp example = ('chr4', '49323694', 'chr22', '12189934')

                #CHECK FIRST GENE
                for fusion in potential_fusions_list: #iterate through the list of possible fusions 
                    # fusion = (gene, mate, chrom, start, end)
                    # temp = (chr, position, mate chr, mate position)
                    if fusion[2] == temp[0] and temp[1] >= fusion[3] and temp[1] <= fusion[4]: #if chromosomes match AND position is between the start and the end of the gene 
                        #CHECK MATE GENE
                        for fusion2 in potential_fusions_list:  #check matched gene
                            if fusion2[0] == fusion[1] and temp[2] == fusion2[2] and temp[3] >= fusion2[3] and temp[3] <= fusion2[4]: #if chromosomes match AND position is between the start and the end of the gene 
                                found_fusion = (patient, fusion[0], fusion2[0], temp) #ex: ('TARGET-20-PAURDN-03A-01D_consensus', 'NUP98', 'NSD1', ('chr11', '3743985', 'chr5', '177233312'))
                                if found_fusion not in found_fusions_list: 
                                    found_fusions_list.append(found_fusion)


#this will print out all the found fusions! 
for fusion in found_fusions_list:
    print(fusion)



#with open(f"{output}/{patient}_CONSENSUS.vcf", "w") as fileout:









#ignore this! 

                    #check second gene (backwards check!)
                    #elif fusion[2] == temp[2] and temp[3] >= fusion[3] and temp[3] <= fusion[4]:
                    #if chromosomes match AND position is between the start and the end of the gene 
                        #matched_gene = fusion[1] 
                        #for fusion2 in potential_fusions_list:  #check matched gene
                            #if fusion2[0] == matched_gene and temp[0] == fusion2[2] and temp[1] >= fusion2[3] and temp[1] <= fusion2[4]:
                            #if chromosomes match AND position is between the start and the end of the gene 
                                #found_fusion = (patient, fusion2[0], fusion[0], temp) #ex: ('TARGET-20-PAURDN-03A-01D_consensus', 'NUP98', 'NSD1', ('chr11', '3743985', 'chr5', '177233312'))
                                #matched_fusion = (patient, fusion2[0], fusion[0], temp) #ex: ('TARGET-20-PAURDN-03A-01D_consensus', 'NSD1', 'NUP98', ('chr11', '3743985', 'chr5', '177233312'))
                                #if found_fusion not in found_fusions_list: #or matched_fusion not in found_fusions_list:
                                    #print("nope!")
                                    #found_fusions_list.append(found_fusion)
