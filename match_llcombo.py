#!/usr/bin/env python
import argparse
import os
import re
import itertools

# ./match_llcombo.py -f /projects/bgmp/shared/groups/2024/aml/lenara/trial_files
  
  
#Functions
def build_infodict(infocol):
    '''Parses through lines in VCf file, pulling and storing the info column (column 7)'''
    splitfo = info.split(";") #splitting info line by ;
    infodict = dict()
    for term in splitfo:
        term = term.split("=")
        #this is items in info like "IMPRECISE" (which have no value) will not get added to the dictionary
        if len(term) > 1:
            infodict[term[0]] = term[1]

    return infodict



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
chrom_list = ("chr1", "chr2", "chr24")
chr_sv_combos =  list(itertools.product(svlist, chrom_list))
# print(chr_sv_combos)


#POPULATING DICTIONARY -----------------------------------------------------------------------------------------------------
#open all files by patient ID 
for patient in patient_list:
    dictionary = {}
    for caller in caller_list:
        dictionary[caller] = {}
    with open(f"{f}/{patient}.{caller_list[0]}.vcf", "r") as file1, open(f"{f}/{patient}.{caller_list[1]}.vcf", "r") as file2, open(f"{f}/{patient}.{caller_list[2]}.vcf", "r") as file3: 
        while True: 
            line1 = file1.readline()
            line_as_list1 = line1.strip().split("\t")
            line2 = file2.readline()
            line_as_list2 = line2.strip().split("\t")
            line3 = file3.readline()
            line_as_list3 = line3.strip().split("\t")

            #caller 1
            if line_as_list1 == ['']:
                break
            
            if line_as_list1[0].startswith("#") != True:
                            
                #get start position
                start_pos = line_as_list1[1]
                #get chr
                chrom = line_as_list1[0]
                #get SVLEN and SVEND
                info = line_as_list1[7]
                infodict = build_infodict(info)
                # print(infodict)
                svtype = infodict['SVTYPE']

                if svtype == "BND":
                #this is where abraham's breakend algorithm will go :)
                    continue
                else:

                    svlen = infodict['SVLEN']
                    end = infodict['END']

                    #get allelic depth
                    #accounts for sniffles having their AD in the info line
                    if "SUPPORT" in infodict:
                        allelic_depth = infodict["SUPPORT"]

                    #if caller is not sniffles
                    else:
                        sample_stats = line_as_list1[9]
                        allelic_depth = sample_stats.split(":")[1]

                    #append to dictionary
                    #if insertion
                    if svtype == "INS":
                        dictionary[caller_list[0]][(chrom, svtype, (start_pos, start_pos+svlen), allelic_depth)] = line1
                    #if not insertion
                    else: 
                        dictionary[caller_list[0]][(chrom, svtype, (start_pos, end), allelic_depth)] = line1
    

            #caller 2
            if line_as_list2 == ['']:
                break 
            
            if line_as_list2[0].startswith("#") != True:
                #get start position
                start_pos = line_as_list2[1]
                #get chr
                chrom = line_as_list2[0]
                #get SVLEN and SVEND
                info = line_as_list2[7]
                infodict = build_infodict(info)
                # print(infodict)
                svtype = infodict['SVTYPE']

                if svtype == "BND":
                #this is where abraham's breakend algorithm will go :)
                    continue
                else:

                    svlen = infodict['SVLEN']
                    end = infodict['END']

                    #get allelic depth
                    #accounts for sniffles having their AD in the info line
                    if "SUPPORT" in infodict:
                        allelic_depth = infodict["SUPPORT"]

                    #if caller is not sniffles
                    else:
                        sample_stats = line_as_list2[9]
                        allelic_depth = sample_stats.split(":")[1]

                    #append to dictionary
                    #if insertion
                    if svtype == "INS":
                        dictionary[caller_list[1]][(chrom, svtype, (start_pos, start_pos+svlen), allelic_depth)] = line2
                    #if not insertion
                    else: 
                        dictionary[caller_list[1]][(chrom, svtype, (start_pos, end), allelic_depth)] = line2

                
            
            #caller 3
            if line_as_list3 == ['']:
                break 
            
            if line_as_list3[0].startswith("#") != True:
                    #get start position
                start_pos = line_as_list3[1]
                #get chr
                chrom = line_as_list3[0]
                #get SVLEN and SVEND
                info = line_as_list3[7]
                infodict = build_infodict(info)
                svtype = infodict['SVTYPE']

                if svtype == "BND":
                #this is where abraham's breakend algorithm will go :)
                    continue
                else:

                    svlen = infodict['SVLEN']
                    end = infodict['END']

                    #get allelic depth
                    #accounts for sniffles having their AD in the info line
                    if "SUPPORT" in infodict:
                        allelic_depth = infodict["SUPPORT"]

                    #if caller is not sniffles
                    else:
                        sample_stats = line_as_list3[9]
                        allelic_depth = sample_stats.split(":")[1]

                    #append to dictionary
                    #if insertion
                    if svtype == "INS":
                        dictionary[caller_list[2]][(chrom, svtype, (start_pos, start_pos+svlen), allelic_depth)] = line3
                    #if not insertion
                    else: 
                        dictionary[caller_list[2]][(chrom, svtype, (start_pos, end), allelic_depth)] = line3
            
        print(dictionary)


