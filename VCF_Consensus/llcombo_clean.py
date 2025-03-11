#!/usr/bin/env python

# Authors: Lena Allen, Lauren Williams, Maura Kautz, Abraham Solomon
# Email: lenarallen19@gmail.com, laurenrwilliams.12@gmail.com, maura.kautz@outlook.com

import argparse
import os
import re
import itertools  
import sys


def overlaps(span1, span2):
    return span1[0] <= span2[1] and span2[0] <= span1[1]

def getoverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def build_infodict(infocol):
    '''Parses through lines in VCFile, pulling and storing the info column (column 7)'''
    splitfo = infocol.split(";") #splitting info line by ;
    infodict = dict()
    for term in splitfo:
        term = term.split("=")
        #this is items in info like "IMPRECISE" (which have no value) will not get added to the dictionary
        if len(term) > 1:
            infodict[term[0]] = term[1]

    return infodict

def BND_parser(line_as_list, infodictionary):
    
    BND_type: int = 0 #Refers to BND type.
    #BND_type = 1: N[
    #BND_type = 2: ]N
    #BND_type = 3: G]
    #BND_type = 4: [C

    #BND_dict = {BND_func_tuple : full line}
    #Tuple of CHROM, POS, MCHROM, MPOS, BND_type

    BND_func_tuple: tuple = ()

    bnd_mate_pos = re.search(r"(?<=:)\d+", line_as_list[4])
    bnd_mate_chrom = re.search(r"\w+(?=:)", line_as_list[4])

    #Extracting AD for breakends.
    if "SUPPORT" in infodictionary:
        allelic_depth = int(infodictionary["SUPPORT"])
    else:
        sample_stats = line_as_list[9]
        #take allelic depth of alternate allele
        allelic_depth = sample_stats.split(":")[1]
        allelic_depth = int(allelic_depth.split(",")[1])

    if re.findall(r"[A,T,G,C,N]+\[", line_as_list[4]):
        # Sample string: G[chr20:29368734[
        # Match: G[
        #MRF
        BND_func_tuple = (line_as_list[0], "BND", 1, allelic_depth, bnd_mate_chrom.group(), bnd_mate_pos.group())
        #print(f"A mate right forward BND call t[p[  {BND_func_tuple}")

    elif re.findall(r"\]+[A,T,G,C,N]", line_as_list[4]):
        # Sample string: ]chr20:29368734]C
        # Match: ]C
        #MLF
        BND_func_tuple = (line_as_list[0], "BND", 2, allelic_depth, bnd_mate_chrom.group(), bnd_mate_pos.group())
        #print(f"A mate left forward BND call ]p]t")

    elif re.findall(r"[A,T,G,C,N]+\]", line_as_list[4]):
        # Sample string: G]chr20:29368734]
        # Match: G]
        #MRR
        BND_func_tuple = (line_as_list[0], "BND", 3, allelic_depth, bnd_mate_chrom.group(), bnd_mate_pos.group())
        #print(f"A mate right reversed BND call t]p]")

    elif re.findall(r"\[+[A,T,G,C,N]", line_as_list[4]):
        # Sample string: [chr20:29368734[C
        # Match: [C
        #MLR
        BND_func_tuple = (line_as_list[0], "BND", 4, allelic_depth, bnd_mate_chrom.group(), bnd_mate_pos.group())
        #print(f"A mate left reversed BND call [p[t")
    return BND_func_tuple
    # BND_func_tuple = (CHROM, SV, BND_type, AD, Mate_CHROM, Mate_POS)
    
def populate_bigdictionary(patientid, callerlist, foldername):
    '''Given a patient ID, a list of callers, and the name of the folder containing the VCF files, 
        returns a dictionary where keys = caller and
        values = {(chrom, svtype, (start_pos, start_pos+svlen), allelic_depth): line}
        for each patient'''
    bigdictionary = {}
    for caller in callerlist:
        bigdictionary[caller] = {}

        with open(f"{foldername}/{patientid}.{caller}.vcf", "r") as file:

            #iterate through lines in VCF file
            while True:
                line = file.readline()
                line_as_list = line.strip().split("\t")

                if line_as_list == ['']:
                    break
                
                #if it's not a header line
                elif line_as_list[0].startswith("#") != True:
                    #get SVTYPE
                    infos = line_as_list[7]
                    infodictionary = build_infodict(infos)
                    svtype = infodictionary['SVTYPE']
                    
                
                    #add caller to each line for later analysis
                    line_as_list[2] = f"{caller}.{svtype}"
                    line = "\t".join(line_as_list)

                    #get start position
                    start_pos = int(line_as_list[1])
                    #get chromosome
                    chrom = line_as_list[0]

                    #get information for breakend lines
                    if svtype == "BND": 
                        mchr = re.search(r"\w+(?=:)", line_as_list[4])
                        if mchr.group() in chrom_list and line_as_list[0] in chrom_list:
                        #parse breakend lines
                            BND_tuple = BND_parser(line_as_list, infodictionary)
                            bigdictionary[caller][BND_tuple] = line

                    else:
                        #get SV length
                        svlen = infodictionary['SVLEN']
                        
                        #if there are two lengths for variants of the same type
                        #with the same start pos, use the longest length
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
                        bigdictionary[caller][(chrom, svtype, (start_pos, start_pos+svlen), allelic_depth)] = line
    return bigdictionary

def filter_duplicates(dict_of_overlaps):
    '''Takes a dictionary where keys = calls from one caller and values = overlapping calls from other callers and 
    1) identifies call with highest allelic depth for each entry
    2) deduplicates entries with same start and end position'''


    #if there are overlaps BETWEEN callers, add the call with the highest allelic depth to highestlist
    #if it's tied, just add the first call
    highestlist = []

    #for each call in the overlap dictionary
    #if the key has the highest allelic depth, append it to highestlist
    #if one of the values (overlapping calls from other callers) has the highest allelic depth,
    #append it to highest list 
    for entry in dict_of_overlaps:
        highestAD = entry
        ad = entry[0][3]

        for i in overlapdict[entry]:
            if i[0][3] > ad:
                highestAD = i
                ad = i[0][3]

        if highestAD not in highestlist:
            highestlist.append(highestAD)
    

    #FILTER DUPLICATES!!! ---------------------------------------------------------
    #initalize set to hold start and end positions (as a tuple) observed for all calls (with the highest AD) found to have overlaps
    #set will only have unique positions
    seen_coordinates = set()
    for highestAD in highestlist:
        seen_coordinates.add(highestAD[0][2])
    
    #initialize empty list to hold de-duplicated list of calls with the highest allelic depths
    topADs = []

    #iterate thorough set of unique positions
    for coord in seen_coordinates:
        #initialize empty list
        coord_list = []
        #iterate through highestlist and add only calls with matching positions
        for highestAD in highestlist:
            if highestAD[0][2] == coord:
                coord_list.append(highestAD)

        #iterate through coord list, appending only calls with highest AD to topADs
        for i, cs in enumerate(coord_list):
            if i == 0:
                topAD = cs
                currentAD = cs[0][3]
            else:
                if cs[0][3] > currentAD:
                    topAD = cs
                    currentAD = cs[0][3]
        topADs.append(topAD)
    return topADs


# ---------------------------------------------------------------------------- #
#                                   ARGPARSE                                   #
# ---------------------------------------------------------------------------- #

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-input", help="Location of folder containing input files", required=True, type=str)
    parser.add_argument("-output", help="Location of folder to place output files", required=True, type=str)
    parser.add_argument("-numcallers", help="Number of callers whose output files are being compared", required=True, type=int)
    parser.add_argument("-min_indels", help="Number of callers that must identify an indel for it to be included in consensus output.", required=True, type=int)
    parser.add_argument("-min_bnds", help="Number of callers that must identify a breakend for it to be included in consensus output.", required=True, type=int)
    parser.add_argument("-min_inv_dup", help="Number of callers that must identify an inversion or duplication for it to be included in consensus output.", required=True, type=int)
    parser.add_argument("-min_overlap", help="Minimum overlap of bases required for two calls from different callers to be considered concordant.", default=1, required=True,type=int)
    return parser.parse_args() 

args = get_args()
input = args.input
output = args.output

numcallers = args.numcallers
min_indels = args.min_indels
min_bnds = args.min_bnds
min_inv_dup = args.min_inv_dup
# keepbnds = args.keep_bnds
minoverlap = args.min_overlap

#other potential options:
#keep decoys?

# ---------------------------------------------------------------------------- #
#                                  PARSE FILES                                 #
# ---------------------------------------------------------------------------- #

#get all file names from folder and save as list
dir_list = os.listdir(input)

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
chrom_list = ("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11","chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21","chr22", "chrX", "chrY", "chrM")
chr_sv_combos =  list(itertools.product(chrom_list, svlist))
sv_mincaller_dict = {"INS": min_indels, "DEL": min_indels, "BND": min_bnds, "INV": min_inv_dup, "DUP": min_inv_dup}

# ---------------------------------------------------------------------------- #
#                             CONSENSUS GENERATION                             #
# ---------------------------------------------------------------------------- #

# --------------------------- DICTIONARY POPULATION -------------------------- #

for patient in patient_list:
    #open output consensus file for current patient
    with open(f"{output}/{patient}_consensus.vcf", "w") as fo:
        
        #populate dictionary where keys = caller, 
        #values = {(chrom, svtype, (start_pos, start_pos+svlen), allelic_depth): line}
        dictionary = populate_bigdictionary(patient, caller_list, input)
        # ---------------------------------FILTERING-------------------------------- #

        #limit search space to one structural variant/chromosome combination
        for combo in chr_sv_combos:
            chrom = combo[0]
            sv = combo[1]
            allbnds = []
                
            #Initialize and populate dictionary to hold only structural variants of currrent structural variant type on current chromosome 
            subsetdict = {}
            for caller in caller_list:
                subsetdict[caller] = {k:v for k,v in dictionary[caller].items() if k[0] == chrom and k[1] == sv}
            
            # ----------------------------- ASSESSING OVERLAP ---------------------------- #
            
            #initialize a dictionary where keys = [(chrom, svtype, (start_pos, start_pos+svlen), allelic_depth), caller] 
            #for every entry in each caller dictionary in subsetdict and values = an empty list
            overlapdict = {}
            for caller in subsetdict:
                for i in subsetdict[caller]:
                    overlapdict[(i, caller)] = []

            

            #populate overlap dictionary
            #for each key, A, in overlapdict, if another key overlaps with A and is not from the same caller, 
            #that key will be appended to A's value, which begins as an empty list
            for entry in overlapdict:
                for other in overlapdict:
                    #if callers are not the same
                    if entry[1] != other[1]:

                        #if structural variant is NOT a breakend
                        if sv != "BND" and overlaps(entry[0][2], other[0][2]) == True:
                            if minoverlap == 1:
                                overlapdict[entry].append(other)
                            else:
                                if getoverlap(entry[0][2], other[0][2]) >= minoverlap:
                                    overlapdict[entry].append(other)

                        #if structural variant IS a breakend
                        elif sv == "BND":
                            if entry[0][4] == other[0][4]:
                                #checking for breakend types
                                if entry[0][2] == other[0][2]:
                                    overlapdict[entry].append(other)
                                elif (entry[0][2] == 2 and other[0][2] == 4) or (entry[0][2] == 3 and other[0][2] == 1):
                                    if entry[0][5] >= other[0][5]:
                                        overlapdict[entry].append(other)
                                elif (entry[0][2] == 4 and other[0][2] == 2) or (entry[0][2] == 1 and other[0][2] == 3):
                                    if entry[0][5] <= other[0][5]:
                                        overlapdict[entry].append(other)
                    
            
            
            #If there are no overlaps between callers, remove that key
            overlapdict = {k:v for k,v in overlapdict.items() if (1 + len(v)) >= sv_mincaller_dict[sv]}

            #filter out duplicates
            tops = filter_duplicates(overlapdict)

            #write to file
            for item in tops:
                line = dictionary[item[1]][item[0]]
                fo.write(f"{line}\n")


                    
                    
                
