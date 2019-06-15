#!/usr/bin/env python
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
    description='Returns a tsv files with detected AMR genes')
parser.add_argument(
    "amr", help="1. Path to file containing abricate results, e.g. AMR.tsv",
    type=argparse.FileType('r'))
parser.add_argument(
    "resf", help="2. Path to file containing resistance info, e.g. resfinder_names.txt",
    type=argparse.FileType('r'))
parser.add_argument(
    "file_out", help="3. Path to the outputfile",
    type=argparse.FileType('w'))
    
args = parser.parse_args()
amr = args.amr
resf = args.resf
file_out = args.file_out

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

################################################################################

# make dictionary out of resfinder names file
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

with file_out as smm:
    # write tsv with counts
    result_writer = csv.writer(smm,delimiter='\t', lineterminator = '\n')
    # header row
    result_writer.writerow(["protein_ID","AMR_gene", "resistance"])
    for k in amr_info.keys():
        # amr_info[k] is a list , iterate over it:
        for gene in amr_info[k]:
            # tet(Q)_1 --> tet
            gene_fam = gene[0:3]
            resfinder.setdefault(gene_fam,"NULL")
            # [0] because dict_maker returns list, list NA in this case, want first element
            ab=resfinder[gene_fam][0]
            line = [k, gene, ab ]
            result_writer.writerow(line)
smm.close()

