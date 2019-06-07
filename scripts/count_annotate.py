#!/exports/cmvm/eddie/eb/groups/watson_grp/software/julie_conda/meta-assembly/bin/python3
# filter blastx results
# column 3 = %Id ; >80
# column 4 = alignment lentgh ; > 30
# create counts file: contig_ID & nr of significant matches

import csv
import sys
import argparse
import re

################################################################################
# VALIDATE INPUT
################################################################################

parser = argparse.ArgumentParser(
    description='Returns several tsv files with kegg info and counts')
parser.add_argument(
    "file_in", help="1. Path to tabular blast results, e.g. matches.m8",
    type=argparse.FileType('r'))
parser.add_argument(
    "info_dir", help="2. Path and prefix to results from kegg_info.py, e.g. annotation/sample/sample_",
    type=str)
parser.add_argument(
    "seq_len", help="3. Sequence length of he raw reads",
    type=int)
    
args = parser.parse_args()
blastx = args.file_in
info = args.info_dir
seq_len = args.seq_len

################################################################################
# create dictionary from file with 2 columns {col1: [col2(a),col2(b)],...}
def dict_maker(file):
    D = {}
    temp=""
    multi=[]
    fileObj = open(file,"r")
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

# combine dictionaries
def dict_combi(dct_counts, dct_in, dct_combi):
    for k,v in dct_counts.items():
        if k in dct_in:
            for _ in dct_in[k]: 
                dct_combi.setdefault(_,0)
                dct_combi[_]+=v

# write a dictionary to a file, sorted descending by value
def dict_writer(dct, filename):
    fileObj=open(filename, mode='w')
    file_writer = csv.writer(fileObj,delimiter='\t', lineterminator = '\n')
    sort_desc=sorted(dct, key=dct.__getitem__, reverse=True)
    for k in sort_desc:
        file_writer.writerow([k,dct[k]])
    fileObj.close()

################################################################################

# set cut-off boundaries:
x = int(seq_len/3) # because nucl to AA seq
lower = x - int(x/10) # for sample with 33 max, chose interval 30-34 ; -10% ? or -3 ?
upper = x + 1 # allow insertion ; again, always one? or depending on length of sequence?

counts={}
print("Counting...")
# read & count results
for hit in blastx.readlines():
    hit = hit.split('\t')
    # col 2 = contig, col3 = %ID , col4 = AA length
    contig=hit[1]
    Id=float(hit[2])
    length=float(hit[3])
    if Id >= 80 and length >= lower and length <= upper :
        # if key contig does not exist yet, value is 0
        counts.setdefault(contig,0)
        counts[contig]+=1
blastx.close()

smm = open(info + "summary.tsv", "r")
# write results to table
print("Combining counts and annotation...")
annotation=info+"counts_summary.tsv"
with open(annotation, mode='w') as result_file:
    result_writer = csv.writer(result_file,delimiter='\t', lineterminator = '\n')
    # header row
    result_writer.writerow(["contig_ID","dmnd_ID", "KO_ID", "name", "EC", "cazy","pathways","counts"])
    for line in smm.readlines():
        line = line.strip().split('\t')
        ref = line[0]
        if ref in counts:
            line.append(counts[ref])
            result_writer.writerow(line)
    smm.close()
result_file.close()

D_enz = dict_maker(info+"enz_acc")
D_cazy = dict_maker(info+"cazy")
D_path = dict_maker(info+"pathways")

C_enz = {}
C_cazy = {}
C_path = {}

dict_combi(counts, D_enz, C_enz)
dict_combi(counts, D_cazy, C_cazy)
dict_combi(counts, D_path, C_path)

dict_writer(counts, info+"counts_protein")
dict_writer(C_enz, info+"counts_enz")
dict_writer(C_cazy, info+"counts_cazy")
dict_writer(C_path, info+"counts_path")

print("Done.")
