#!/usr/bin/env python
import argparse
import os
import re
import itertools

# ./canonfusions.py -input /projects/bgmp/shared/groups/2024/aml/shared/llcombo10_output_files -output /projects/bgmp/shared/groups/2024/aml/lwil/canon_fusions
  
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
chrom_list = ("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11","chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21","chr22")
chr_sv_combos =  list(itertools.product(svlist, chrom_list))

BND_dict: dict = {}
SV1_BND_dict: dict = {}
SV2_BND_dict: dict = {}

#open all files by patient ID 
for patient in patient_list:
    dictionary = {}
    for caller in caller_list:
        dictionary[caller] = {}
    #POPULATING DICTIONARY -----------------------------------------------------------------------------------------------------
    with open(f"{input}/{patient}.{caller_list[0]}.vcf", "r") as file1, open(f"{input}/{patient}.{caller_list[1]}.vcf", "r") as file2, open(f"{input}/{patient}.{caller_list[2]}.vcf", "r") as file3:
        while True:
            line1 = file1.readline()
            line_as_list1 = line1.strip().split("\t")
            line2 = file2.readline()
            line_as_list2 = line2.strip().split("\t")
            line3 = file3.readline()
            line_as_list3 = line3.strip().split("\t")

            #caller 1-----------------------------------------------------------------------------------------------------------
            if line_as_list1 == ['']:
                break
            
            elif line_as_list1[0].startswith("#") != True:
                line_as_list1[2] = caller_list[0]
                #get start position
                start_pos = int(line_as_list1[1])
                #get chr
                chrom = line_as_list1[0]
                #get SVLEN and SVEND
                info = line_as_list1[7]
                infodict = build_infodict(info)
                svtype = infodict['SVTYPE']

                if svtype == "BND":
                    #this is where abraham's breakend algorithm will go :)
                    BND_tuple = BND_parser(line_as_list1)
                    print(BND_tuple)
                    SV1_BND_dict[BND_tuple] = line1
                    continue
                else:

                    svlen = infodict['SVLEN']
                    if "," in svlen:
                        svlen = svlen.split(",")
                    else:
                        svlen = abs(int(svlen))

                    #get allelic depth
                    #accounts for sniffles having their AD in the info line
                    if "SUPPORT" in infodict:
                        allelic_depth = infodict["SUPPORT"]

                    #if caller is not sniffles
                    else:
                        sample_stats = line_as_list1[9]
                        allelic_depth = sample_stats.split(":")[1]
                        allelic_depth = allelic_depth.split(",")[1]

                    #append to dictionary
                    if type(svlen) != list:
                        dictionary[caller_list[0]][(chrom, svtype, (start_pos, start_pos+svlen), allelic_depth)] = line1
                    else:
                        svlen = [int(x) for x in svlen]
                        maxlen = max(svlen)
                        dictionary[caller_list[0]][(chrom, svtype, (start_pos, start_pos+maxlen), allelic_depth)] = line1

            #caller 2-----------------------------------------------------------------------------------------------------------
            if line_as_list2 == ['']:
                break 
            
            elif line_as_list2[0].startswith("#") != True:
                line_as_list2[2] = caller_list[1]

                #get start position
                start_pos = int(line_as_list2[1])
                #get chr
                chrom = line_as_list2[0]
                #get SVLEN and SVEND
                info = line_as_list2[7]
                infodict = build_infodict(info)
                svtype = infodict['SVTYPE']

                if svtype == "BND":
                    #this is where abraham's breakend algorithm will go :)
                    BND_tuple = BND_parser(line_as_list1)
                    print(BND_tuple)
                    SV1_BND_dict[BND_tuple] = line1
                    continue
                else:

                    svlen = infodict['SVLEN']
                    if "," in svlen:
                        svlen = svlen.split(",")
                    else:
                        svlen = abs(int(svlen))
                    end = infodict['END']

                    #get allelic depth
                    #accounts for sniffles having their AD in the info line
                    if "SUPPORT" in infodict:
                        allelic_depth = infodict["SUPPORT"]

                    #if caller is not sniffles
                    else:
                        sample_stats = line_as_list2[9]
                        allelic_depth = sample_stats.split(":")[1]
                        allelic_depth = allelic_depth.split(",")[1]

                    #append to dictionary
                    if type(svlen) != list:
                        dictionary[caller_list[1]][(chrom, svtype, (start_pos, start_pos+svlen), allelic_depth)] = line2
                    else:
                        svlen = [int(x) for x in svlen]
                        maxlen = max(svlen)
                        dictionary[caller_list[1]][(chrom, svtype, (start_pos, start_pos+maxlen), allelic_depth)] = line2


            #caller 3-----------------------------------------------------------------------------------------------------------
            if line_as_list3 == ['']:
                break 
            
            elif line_as_list3[0].startswith("#") != True:
                line_as_list3[2] = caller_list[2]

                #get start position
                start_pos = int(line_as_list3[1])
                #get chr
                chrom = line_as_list3[0]
                #get SVLEN and SVEND
                info = line_as_list3[7]
                infodict = build_infodict(info)
                svtype = infodict['SVTYPE']

                if svtype == "BND":
                    #this is where abraham's breakend algorithm will go :)
                    BND_tuple = BND_parser(line_as_list1)
                    print(BND_tuple)
                    SV1_BND_dict[BND_tuple] = line1
                    continue
                else:

                    svlen = infodict['SVLEN']
                    if "," in svlen:
                        svlen = svlen.split(",")
                    else:
                        svlen = abs(int(svlen))
                    end = infodict['END']

                    #get allelic depth
                    #accounts for sniffles having their AD in the info line
                    if "SUPPORT" in infodict:
                        allelic_depth = infodict["SUPPORT"]

                    #if caller is not sniffles
                    else:
                        sample_stats = line_as_list3[9]
                        allelic_depth = sample_stats.split(":")[1]
                        allelic_depth = allelic_depth.split(",")[1]

                    #append to dictionary
                    if type(svlen) != list:
                        dictionary[caller_list[2]][(chrom, svtype, (start_pos, start_pos+svlen), allelic_depth)] = line3
                    else:
                        svlen = [int(x) for x in svlen]
                        maxlen = max(svlen)
                        dictionary[caller_list[2]][(chrom, svtype, (start_pos, start_pos+maxlen), allelic_depth)] = line3

    with open(f"{output}/{patient}_CONSENSUS.vcf", "w") as fo:
        for callers in caller_list:
            for key in dictionary[caller]:
                fo.write(dictionary[caller][key])