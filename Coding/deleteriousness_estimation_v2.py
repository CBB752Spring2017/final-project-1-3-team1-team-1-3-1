#!/usr/bin/python

__author__ = "Zhaolong Yu"
__copyright__ = "Copyright 2017"
__credits__ = ["Zhaolong Yu"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Zhaolong Yu"
__email__ = "zhaolong.yu@yale.edu"

### Usage:      python3 deleteriouness_estimation_v2.py -i <input file>
### Example:    python3 deleteriouness_estimation_v2.py -i input.txt


import argparse
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import networkx as nx
import mygene
import seaborn

## open the file
parser = argparse.ArgumentParser(description='Centrality Calculation V1.0')
parser.add_argument('-i', '--input', help='input mutation file', required=True)

args = parser.parse_args()
#os.chdir("/Users/jerome/Projects/deleteriousness_estimation")
fo_snp = open(args.input,"r")

## build the snp info dict
snv_dict = {}
lc = 0
fo_snp.seek(0)
for rawline in fo_snp:
    line = rawline.strip().split("\t")
    
    ## extract the gene
    gene = line[2].split(";")[1]
    
    ## generate the id
    snp_id = line[4] + "_" + line[1]
    
    ## mutation type
    mut_status = line[3]
    
    ## mutation details
    snp_info = line[1].split(":")[2:4]
    aa_info = line[5].split("_")[3].split("->")
    if len(aa_info[1]) > 1:
        aa_info[1] = aa_info[1][0]
    
    ## gerp score
    gerp = line[2].split(";")[0].split("=")[1]
    
    if snp_id not in snv_dict.keys():
        snv_dict[snp_id] = []
        snv_dict[snp_id] = [[line[1],gene,mut_status,gerp,snp_info,aa_info,line[6],line[7]]]
    else:
        snv_dict[snp_id].append([[line[1],gene,mut_status,gerp,snp_info,aa_info,line[6],line[7]]])
    lc = lc + 1
fo_snp.close()

id_list = []

for item in snv_dict.keys():
    id_list.append(item)
    
## the transcript in total
len(snv_dict.keys())

count_multi = 0
multi_list = []
for transcript in snv_dict.keys():
    if len(snv_dict[transcript]) > 1:
        count_multi = count_multi + 1


## summarize different kinds of mutations
fo_snp = open(args.input,"r")
fo_snp.seek(0)
syn = []
non_syn = []
premature = []
indels = []
others = []
lc = 0
for rawline in fo_snp:
    line = rawline.strip().split("\t")
    syn_status = line[3]
    snp_info = line[1].split(":")[2:4]
    if syn_status == "nonsynonymous":
        aa_info = line[5].split("_")[3].split("->")
        non_syn.append([line[4],snp_info,aa_info])        
    elif syn_status == "synonymous":
        syn.append([line[4],snp_info])
    elif syn_status == "prematureStop":
        ori_pep = line[6]
        new_pep = line[7].split("*")[0]
        length_dif = len(ori_pep)-len(new_pep)
        premature.append([line[4],snp_info,length_dif])
    elif syn_status == "indels":
        if snp_info[1] == "*":
            if len(snp_info[0])%3 != 0:
                indels.append([line[4],5.5,"del_frameshift"])
            else:
                indels.append(line[4],len(snp_info[0])/3, "deletion")
        else:
            if (len(snp_info[1])-len(snp_info[0]))%3 != 0:
                indels.append([line[4],5.5,"insert_frameshift"])
            else:
                indels.append(line[4],(len(snp_info[1])-len(snp_info[0]))/3, "insertion")
    else:
        others.append([line[3],line[4]])
        print(line)
    lc = lc + 1
fo_snp.close()

## deleteriousness estimation
## read the files

fo_pam250 = open("pam250.txt","r")
fo_codon_usage = open("codon_usage.txt","r")
fo_aa = open('scaled_aa_info2.txt',"r")

## build a aa_info dict including all the scaled properties of 20 amino acids
fo_aa.readline()
aa_dict = {}
for line in fo_aa:
    new_line = line.strip().split()
    vector = []
    for item in new_line[1:]:
        vector.append(float(item))    
    aa_dict[new_line[0]] = np.array(vector)
    

fo_aa.close()

## build a dict to store the usage bias for different codons
fo_codon_usage.seek(0)
codon_dict = {}
for rawline in fo_codon_usage:
    line = rawline.strip().split(")")
    for short_str in line[:4]:
        subline = short_str.strip().split(" ")
        if subline[1] not in codon_dict.keys():
            codon_dict[subline[1]] = []
            codon_dict[subline[1]] = [[subline[0],subline[0][2],subline[2]]]
        else:
            codon_dict[subline[1]].append([subline[0],subline[0][2],subline[2]])
 
## build a dict to store the PAM250 matrix           
fo_pam250.seek(0)
fo_pam250.readline()    
pam250_dict = {}        
aa = fo_pam250.readline().strip().split()
for i in aa:
    pam250_dict[i] = {}
for i in range(20):
    line = fo_pam250.readline().strip().split()[1:]
    sub_aa = line[0]
    vector = line[1:]
    for j in range(20):
        pam250_dict[aa[j]][sub_aa] = float(vector[j])/100

vector_list = []        
np_id = 0
        
def deleteriousness_estimation(snv_element):
    
    ## extract the info
    mut_status = snv_element[0][2]
    snp_info = snv_element[0][4]
    aa_info = snv_element[0][5]
    line = snv_element[0]
    tmp_id = ":".join(line[0].split(":")[0:2])
    snps = "".join(line[0].split(":")[2:4])
    tmp_id = tmp_id + "-" + snps
    gerp = 0
    if snv_element[0][3] != "." :
        if float(snv_element[0][3]) > 6: 
            gerp = float(snv_element[0][3])
    if mut_status == "nonsynonymous":
        PAM_prob = pam250_dict[snv_element[0][5][0]]
        standard_vector = aa_dict[snv_element[0][5][0]]
        substitute_vector = aa_dict[snv_element[0][5][1]]
        resi_vector = abs(substitute_vector - standard_vector)
        deviance_snv =  sum(abs(substitute_vector - standard_vector))/150
        est_score = 0.05 + deviance_snv + gerp/10     
    elif mut_status == "synonymous":
        aa = aa_info[0]
        original_codon = snp_info[0]
        new_codon = snp_info[1]
        usage = codon_dict[aa]
        for codon in usage:
            if original_codon == codon[0]:
                orginal_usage = codon[2]
            if new_codon ==codon[0]:
                orginal_usage = codon[2]
        if original_usage < new_usage:
            est_score = 0
        else:
            est_score = 0.01 * (original_usage - new_usage)/original_usage
      
    elif mut_status == "prematureStop":
        ori_pep = snv_element[0][6]
        new_pep = snv_element[0][7].split("*")[0]
        length_dif = len(ori_pep)-len(new_pep)
        if length_dif <= 5:
            est_score = 0.8 + 0.4 * length_dif
        else:
            est_score = 1   
    elif mut_status == "indels":
        if snp_info[1] == "*":
            # frameshift caused by deletion
            if len(snv_element[0][4][0])%3 != 0:
                est_score = 1
            # no frameshift
            else:
                est_score = 0.8 + 0.4 * len(snv_element[0][4][0])/3
        else:
            # frameshift caused by insertion
            if (len(snp_info[1])-len(snp_info[0]))%3 != 0:
                est_score = 1
            # no frameshift
            else:
                est_score = 0.8 + 0.4 * len(snv_element[0][4][0])/3
    else:
        est_score = 0
    return(est_score)

## build the score dict
deleteriousness_est = {}
score_est = {}
## calculate the deleteriousness score
index = 0
for snv in id_list:
    tmp_score = deleteriousness_estimation(snv_dict[snv])
    deleteriousness_est[snv] = [snv_dict[snv][0][1],snv_dict[snv][0][2],tmp_score]
    score_est[snv] = tmp_score
    index = index +1 

## sort the score
   
sorted_keys = sorted(score_est, key=score_est.get, reverse=True)

## output the deleteriousness score    
fo_output = open("deleterious_score.txt","w")
for snv in sorted_keys:
    fo_output.write(snv + "\t" + "\t".join(deleteriousness_est[snv][:2]) + "\t" + str(deleteriousness_est[snv][2]) + "\n")
fo_output.close() 

