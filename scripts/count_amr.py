#!/exports/cmvm/eddie/eb/groups/watson_grp/software/julie_conda/meta-assembly/bin/python3
# arg 1) file with proteincounts "kxxx   xxxx" 
# arg 2) abricate results

import csv
import sys
import argparse
import re

################################################################################
# VALIDATE INPUT
################################################################################

parser = argparse.ArgumentParser(
    description='Returns a tsv files with kegg info and name')
parser.add_argument(
    "counts", help="1. Path to protein counts, e.g. counts_protein.tsv",
    type=argparse.FileType('r'))
parser.add_argument(
    "amr", help="2. Path to file containing abricate results, e.g. AMR.tsv",
    type=argparse.FileType('r'))
parser.add_argument(
    "resf", help="3. Path to file containing resistance info, e.g. resfinder_names.txt",
    type=argparse.FileType('r'))
parser.add_argument(
    "file_out", help="4. Outputfile",
    type=argparse.FileType("w"))
    
args = parser.parse_args()
cts = args.counts
amr = args.amr
resf = args.resf
file_out = args.file_out
#"C:/Users/- Julie -/AppData/Local/Programs/Python/Python37/python.exe" "c:/Users/- Julie -/Documents/BIT/Roslin stage/project/amr.py"
#cts = open("c:/Users/- Julie -/Documents/BIT/Roslin stage/project/10_annotation/chicken_ERR2241690/chicken_ERR2241690_counts_protein", "r")
#amr = open("c:/Users/- Julie -/Documents/BIT/Roslin stage/project/chicken_ERR2241690_AMR.tsv", "r")
#resf = open("c:/Users/- Julie -/Documents/BIT/Roslin stage/project/resfinder_names.txt", "r")
#file_out = open("c:/Users/- Julie -/Documents/BIT/Roslin stage/project/chicken_ERR2241690_AMR_counts.tsv", "w")
################################################################################
# Functions
################################################################################

# create dictionary from file with 2 columns {col1: [col2(a),col2(b)],...}
def dict_maker(fileObj):
    D = {}
    temp=""
    multi=[]
    lines = fileObj.readlines()
    lines.sort()
    for line in lines:
        line=line.split('\t')
        if temp != line[0]:
            multi=[line[1].strip()]
            D[line[0]] = multi
        else:
            multi.append(line[1].strip())
            D[line[0]] = multi
        temp=line[0]
    fileObj.close()
    return D

# write a dictionary to a file, sorted descending by value
def dict_writer(dct, filename):
    #fileObj=open(filename, mode='w')
    # file already opened by argparse
    fileObj = filename
    file_writer = csv.writer(fileObj,delimiter='\t', lineterminator = '\n')
    sort_desc=sorted(dct, key=dct.__getitem__, reverse=True)
    for k in sort_desc:
        file_writer.writerow([k,dct[k]])
    fileObj.close()

# dct x with info about kx_xxx , dict y with counts for kx_xxx_xx
# sum of counts for parts of protein
def dict_count_combi(x,y):
    counts = {}
    for i in x.keys():
        for j in y.keys():
            if re.search( (i +"_") ,j):
                counts.setdefault(i,0)
                counts[i]+= int(y[j][0])
    return counts

################################################################################

# make dictionary out of protein counts:
protein_counts = dict_maker(cts)
resfinder = dict_maker(resf)

# make dictionary with AMR info, storing { 'kxxx': ['gene', 'gene', ..],...}
# edit dict_maker for abricate results
amr_info = {}
temp=""
multi=[]
amr_list = amr.readlines()
# remove first line with headers:
amr_list.pop(0)
for line in amr_list:
    line=line.split('\t')
    seq=line[1]
    gene=line[4]
    cov=float(line[8])
    Id=float(line[9])
    if Id >= 80 and cov >= 80:
        if temp != seq:
            multi=[gene.strip()]
            amr_info[seq] = multi
        else:
            multi.append(gene.strip())
            amr_info[seq] = multi
    temp=seq

# sum counts for parts of protein
counts = dict_count_combi(amr_info,protein_counts)

# write tsv with counts
result_writer = csv.writer(file_out,delimiter='\t', lineterminator = '\n')
# header row
result_writer.writerow(["protein_ID","AMR_gene", "resistance", "abundance"])
for k in amr_info.keys():
    for gene in amr_info[k]:
        gene_fam = gene.split('_')[0]
        line = [k, gene,resfinder[gene_fam][0], counts[k] ]
        result_writer.writerow(line)
file_out.close()
