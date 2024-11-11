#!/usr/bin/env python
import argparse
import os
import re
import itertools

# ./llcombo4_no_filter_within_caller.py -f /projects/bgmp/shared/groups/2024/aml/lenara/shortfiles  
  
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


def BND_parser(line_as_list, infodictionary):
    #NOV 10 --- NEED TO IMPLEMENT THIS:
        #If chrom is in chrom list (a real chromosome) and if mate chrom is in chrom list (a real chromosome) to remove unknown chromosomes.
    
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

    #print(bnd_mate_pos.group())

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

    # BND_func_tuple = (CHROM, SV, BND_type, AD, Mate_CHROM, Mate_POS)
    return BND_func_tuple

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
                        BND_tuple = BND_parser(line_as_list, infodictionary)
                        bigdict[caller][BND_tuple] = line

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
    parser.add_argument("-i", help="Location of folder containing input files", required=True, type=str)
    parser.add_argument("-o", help="Location of folder to place output files", required=True, type=str)
    return parser.parse_args() 

args = get_args()
i = args.i
o = args.o
#get all file names from folder and save as list
dir_list = os.listdir(i)

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
chrom_list = ("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11","chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21","chr22")
chr_sv_combos =  list(itertools.product(chrom_list, svlist))

#open all files by patient ID 
for patient in patient_list:
    with open(f"{o}/testconsensus_no_filter_deduped.vcf", "w") as fo:
        
        #populate dictionary where keys = caller, 
        #values = {(chrom, svtype, (start_pos, start_pos+svlen), allelic_depth): line}
        dictionary = populate_bigdict(patient, caller_list, i)

        #Filtering and comparison --- a lot of this can be turned into a function
        #Right now, this will only handle 3 callers - going to try to fix this over the weekend
        
        #filter down to one sv chrom combo
        for combo in chr_sv_combos:
            chrom = combo[0]
            sv = combo[1]
            subsetdict = {}
            for caller in caller_list:
            #initialize two dictionaries
                subsetdict[caller] = {k:v for k,v in dictionary[caller].items() if k[0] == chrom and k[1] == sv}
            
            overlapdict = {}
            for caller in subsetdict:
                for i in subsetdict[caller]:
                    overlapdict[(i, caller)] = []

            for entry in overlapdict:
                for other in overlapdict:
                    if entry[1] != other[1]:
                        if sv != "BND" and overlaps(entry[0][2], other[0][2]) == True:
                            overlapdict[entry].append(other)
                        elif sv == "BND" and (entry[0][4] == other[0][4]):
                            if entry[0][2] == other[0][2]:
                                overlapdict[entry].append(other)
                            elif (entry[0][2] == 2 and other[0][2] == 4) or (entry[0][2] == 3 and other[0][2] == 1):
                                if entry[0][5] >= other[0][5]:
                                    overlapdict[entry].append(other)
                            elif (entry[0][2] == 4 and other[0][2] == 2) or (entry[0][2] == 1 and other[0][2] == 3):
                                if entry[0][5] <= other[0][5]:
                                    overlapdict[entry].append(other)
            
            #If there are no overlaps between callers, remove that key
            overlapdict = {k:v for k,v in overlapdict.items() if v}

            #if there are overlaps between callers, write the call with the highest AD to file
            #if it's a tie, just print the first one
            highestlist = []
            for entry in overlapdict:
                highestAD = entry
                ad = entry[0][3]
                for i in overlapdict[entry]:
                    if i[0][3] > ad:
                        highestAD = i
                if highestAD not in highestlist:
                    highestlist.append(highestAD)


            for highestAD in highestlist:
                print(f"{entry}\t{highestAD[0][2]}\t{highestAD[1]}")
                line = dictionary[highestAD[1]][highestAD[0]]
                fo.write(f"{line}\n")
            
