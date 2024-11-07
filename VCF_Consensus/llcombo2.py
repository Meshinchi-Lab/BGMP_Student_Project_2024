#!/usr/bin/env python
import argparse
import os
import re
import itertools

# ./llcombo2.py -f /projects/bgmp/shared/groups/2024/aml/lenara/trial_files
  
#Functions
def overlaps(span1, span2):
    return span1[0] <= span2[1] and span2[0] <= span1[1]

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
    parser.add_argument("-input", help="Location of folder containing input files", required=True, type=str)
    parser.add_argument("-output", help="Location of folder where the output file will go", required=True, type=str)
    return parser.parse_args() 

args = get_args()
input = args.input
output = args.output

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
chrom_list = ("chr1", "chr2")
chr_sv_combos =  list(itertools.product(svlist, chrom_list))

#open all files by patient ID 
for patient in patient_list:
    dictionary = {}
    for caller in caller_list:
        dictionary[caller] = {}
    #POPULATING DICTIONARY -----------------------------------------------------------------------------------------------------
    with open(f"{input}/{patient}.{caller_list[0]}.vcf", "r") as file1, open(f"{input}/{patient}.{caller_list[1]}.vcf", "r") as file2, open(f"{input}/{patient}.{caller_list[2]}.vcf", "r") as file3, open(f"{output}/{patient}.consensus.vcf", "w") as fo:
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

        #FILTERING -----------------------------------------------------------------------------------------------------
        for chromosome in chrom_list:
            for sv in svlist:
                if sv in ["INS", "DEL"]:
                    combodict1 = {}
                    combodict2 = {}
                    for caller in dictionary:
                        combodict1[caller] = {k:v for k,v in dictionary[caller].items() if k[0] == chromosome and k[1] == sv}
                        combodict2[caller] = {k:v for k,v in dictionary[caller].items() if k[0] == chromosome and k[1] == sv}
                        
                        #filter out overlaps within callers
                        #i.e. if two insertions from the same caller overlap
                        #keep the one with the higher AD and remove the one
                        #with lower AD from combodict2
                        for current_key in combodict1[caller]:
                            for key in combodict1[caller]:
                                if key != current_key:
                                    if overlaps(current_key[2], key[2]) == True:
                                        ck_AD = current_key[3]
                                        k_AD = key[3]
                                        if ck_AD >= k_AD and key in combodict2[caller]:
                                            del combodict2[caller][key]
                                        elif k_AD >= ck_AD and current_key in combodict2[caller]:
                                            del combodict2[caller][current_key]
                        
                        #COMPARING ----------------------------------------------------------------------------------------------------
                        
                        if len(combodict2) == 2:
                            lx = list(combodict2.keys())[0]
                            ly = list(combodict2.keys())[1]
                            for tupx in combodict2[lx]:
                                spanx = tupx[2]
                                adx = tupx[3]
                                for tupy in combodict2[ly]:
                                    spany = tupy[2]
                                    ady = tupy[3]
                                    if overlaps(spanx, spany) == True:
                                        if adx >= ady:
                                            fo.write(f"{combodict2[lx][tupx]}\n")
                                        else:
                                            fo.write(f"{combodict2[ly][tupy]}\n")
                                        break

                        elif len(combodict2) == 3:
                            l1 = list(combodict2[caller_list[0]].keys())
                            l2 = list(combodict2[caller_list[1]].keys())
                            l3 = list(combodict2[caller_list[2]].keys())
                            for tup1 in l1:
                                span1 = tup1[2]
                                ad1 = tup1[3]
                                for tup2 in l2:
                                    span2 = tup2[2]
                                    ad2 = tup2[3]
                                    for tup3 in l3:
                                        span3 = tup3[2]
                                        ad3 = tup3[3]
                                        winnerdict = {k:None for k in caller_list}
                                        # is there overlap?
                                        if overlaps(span1, span2) == True:
                                            #which has highest AD?
                                            if ad1 >= ad2:
                                                winnerdict[caller_list[0]] = tup1
                                            else:
                                                winnerdict[caller_list[1]] = tup2
                                        # is there overlap?
                                        if overlaps(span1, span3) == True:
                                            if ad1 >= ad3:
                                                winnerdict[caller_list[0]] = tup1
                                            else:
                                                winnerdict[caller_list[2]] = tup3
                                        
                                        if overlaps(span2, span3) == True:
                                            if ad2 >= ad3:
                                                winnerdict[caller_list[1]] = tup2
                                            else:
                                                winnerdict[caller_list[2]] = tup3
                                        
                                        #remove keys with empty values
                                        winnerdict = {k: v for k,v in winnerdict.items() if v}

                                        if len(winnerdict) == 2:
                                            for currentcaller in winnerdict:
                                                cad1 = winnerdict[currentcaller][3]
                                                for othercaller in winnerdict:
                                                    if currentcaller != othercaller:
                                                        kad1 = winnerdict[othercaller][3]
                                                        if cad1 >= kad1:
                                                            fo.write(f"{combodict2[currentcaller][winnerdict[currentcaller]]}\n")
                                                        else:
                                                            fo.write(f"{combodict2[othercaller][winnerdict[othercaller]]}\n")
                                                        break
                                                break
                                        
                                        elif len(winnerdict) == 1:
                                            for currentcaller in winnerdict:
                                                fo.write(f"{combodict2[currentcaller][winnerdict[currentcaller]]}\n")
                                                

                if sv in ["INV", "DUP"]:
                    # print("invdup")
                    combodict1 = {}
                    combodict2 = {}
                    for caller in dictionary:
                        if caller != caller_list[0]:
                            combodict1[caller] = {k:v for k,v in dictionary[caller].items() if k[0] == chromosome and k[1] == sv}
                            combodict2[caller] = {k:v for k,v in dictionary[caller].items() if k[0] == chromosome and k[1] == sv}

                            for current_key in combodict1[caller]:
                                for key in combodict1[caller]:
                                    if key != current_key:
                                        if overlaps(current_key[2], key[2]) == True:
                                            ck_AD = current_key[3]
                                            k_AD = key[3]
                                            if ck_AD >= k_AD and key in combodict2[caller]:
                                                del combodict2[caller][key]
                                            elif k_AD >= ck_AD and current_key in combodict2[caller]:
                                                del combodict2[caller][current_key]

                            
                                combodict2 = {k: v for k,v in combodict2.items() if len(v)!= 0}
                                
                                if len(combodict2) > 1:
                                    l2 = list(combodict2[caller_list[1]].keys())
                                    l3 = list(combodict2[caller_list[2]].keys())
                                    
                                    for tup2 in l2:
                                        span2 = tup2[2]
                                        ad2 = tup2[3]
                                        for tup3 in l3:
                                            span3 = tup3[2]
                                            ad3 = tup3[3]
                                            if overlaps(span2, span3) == True:
                                                if ad2 >= ad3:
                                                    fo.write(f"{combodict2[caller_list[1]][tup2]}\n")
                                                else:
                                                    fo.write(f"{combodict2[caller_list[2]][tup3]}\n")
                                                break
