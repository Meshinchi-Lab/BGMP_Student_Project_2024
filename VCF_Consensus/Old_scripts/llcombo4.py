#!/usr/bin/env python
import argparse
import os
import re
import itertools

# ./llcombo4.py -f /projects/bgmp/shared/groups/2024/aml/lenara/shortfiles  
  
#Functions
def overlaps(span1, span2):
    return span1[0] <= span2[1] and span2[0] <= span1[1]

def build_infodict(infocol):
    '''Parses through lines in VCf file, pulling and storing the info column (column 7)'''
    splitfo = infocol.split(";") #splitting info line by ;
    infodict = dict()
    for term in splitfo:
        term = term.split("=")
        #this is items in info like "IMPRECISE" (which have no value) will not get added to the dictionary
        if len(term) > 1:
            infodict[term[0]] = term[1]

    return infodict

def populate_bigdict(patientid, callerlist, foldername):
    bigdict = {}
    for caller in callerlist:
        bigdict[caller] = {}
        with open(f"{foldername}/{patientid}.{caller}.vcf", "r") as file:
            while True:
                line = file.readline()
                line_as_list = line.strip().split("\t")

                if line_as_list == ['']:
                    break

                elif line_as_list[0].startswith("#") != True:
                    #add caller to each line for later analysis
                    line_as_list[2] = caller
                    line = "\t".join(line_as_list)

                    #get start position
                    start_pos = int(line_as_list[1])
                    #get chr
                    chrom = line_as_list[0]

                    #get SVTYPE
                    infos = line_as_list[7]
                    infodictionary = build_infodict(infos)
                    svtype = infodictionary['SVTYPE']

                    #temporary diversion for breakends
                    if svtype == "BND":
                        #this is where abraham's breakend algorithm will go :)
                        continue
                    else:
                        #get SV length
                        svlen = infodictionary['SVLEN']
                        
                        #if there are two lengths for variants of the same type
                        #with the same start pos
                        if "," in svlen:
                            svlen = svlen.split(",")
                            svlen = [abs(int(x)) for x in svlen]
                            svlen = max(svlen)
                        else:
                            svlen = abs(int(svlen))
                        
                        #get allelic depth
                        #accounts for sniffles having their AD in the info line
                        if "SUPPORT" in infodictionary:
                            allelic_depth = int(infodictionary["SUPPORT"])

                        #if caller is not sniffles
                        else:
                            sample_stats = line_as_list[9]
                            #take allelic depth of alternate allele
                            allelic_depth = sample_stats.split(":")[1]
                            allelic_depth = int(allelic_depth.split(",")[1])
                        
                        #append to dictionary
                        bigdict[caller][(chrom, svtype, (start_pos, start_pos+svlen), allelic_depth)] = line
    return bigdict

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

#open all files by patient ID 
for patient in patient_list:
    with open(f"testconsensus4.vcf", "w") as fo:
        
        #populate dictionary where keys = caller, 
        #values = {(chrom, svtype, (start_pos, start_pos+svlen), allelic_depth): line}
        dictionary = populate_bigdict(patient, caller_list, f)

        #Filtering and comparison --- a lot of this can be turned into a function
        #Right now, this will only handle 3 callers - going to try to fix this over the weekend
        
        #filter down to one sv chrom combo
        for combo in chr_sv_combos:
            sv = combo[0]
            chrom = combo[1]

            #initialize two dictionaries
            combodict = {}
            combodict2 = {caller:{} for caller in caller_list}
            for caller in caller_list:

                
                #Populate c_overlapdict
                #A dictionary to store calls that overlap WITHIN callers
                #keys = each call, values = [calls that overlap within the same caller]
                c_overlapdict = {}
                #combodict = dictionary with all lines from main dictionary from the same
                #chromosome and of the same sv type
                combodict[caller] = {k:v for k,v in dictionary[caller].items() if k[0] == chrom and k[1] == sv}
                for key in combodict[caller]:
                    c_overlapdict[key] = []
                    for otherkey in combodict[caller]:
                        if otherkey != key:
                            if overlaps(key[2], otherkey[2]) == True:
                                c_overlapdict[key].append(otherkey)
                
                #Filter out overlaps within callers:
                for key in c_overlapdict:
                    #if no overlaps within caller, add call to combodict2
                    if len(c_overlapdict[key]) == 0:
                        combodict2[caller][key] = []
                    #if there are overlaps within caller
                    elif len(c_overlapdict[key]) > 0:
                        #only add the one with the highest AD to combodict2
                        highest = key
                        ad = key[3]
                        for x in c_overlapdict[key]:
                            xad = x[3]
                            if xad > ad:
                                highest = x
                        if highest not in combodict2[caller]:
                            combodict2[caller][highest] = []
            
            #initialize and populate overlapdict
            #A dictionary to store overlaps BETWEEN callers
            #keys = each call, values = [calls that overlap from different callers]
            overlapdict = {}
            for caller in combodict2:
                for i in combodict2[caller]:
                    overlapdict[(i, caller)] = []

            for entry in overlapdict:
                for other in overlapdict:
                    if entry[1] != other[1]:
                        if overlaps(entry[0][2], other[0][2]) == True:
                            overlapdict[entry].append(other)
            
            #If there are no overlaps between callers, remove that key
            overlapdict = {k:v for k,v in overlapdict.items() if v}
            

            #if there are overlaps between callers, write the call with the highest AD to file
            #if it's a tie, just print the first one
            for entry in overlapdict:
                highest = entry
                ad = entry[0][3]
                for i in overlapdict[entry]:
                    if i[0][3] > ad:
                        highest = i
                print(f"{highest[0][2]}\t{highest[1]}")
                line = dictionary[highest[1]][highest[0]]
                fo.write(f"{line}\n")