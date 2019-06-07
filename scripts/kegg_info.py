#!/exports/cmvm/eddie/eb/groups/watson_grp/software/julie_conda/meta-assembly/bin/python3

import csv
import re
import argparse
import os
import sys
#from pymongo import MongoClient

################################################################################
# VALIDATE INPUT
################################################################################

parser = argparse.ArgumentParser(
    description='Create several files containing annotation and a summary.')
parser.add_argument(
    "file_in", help="1. Path to tabular blast results, e.g. diamond/matches.m8",
    type=argparse.FileType('r'))
parser.add_argument(
    "path_annot", help="2. Path to directory containing kegg info flat files",
    type = PathType(exists=True, type='dir'))
parser.add_argument(
    "dir_out", help="3. Prefix for the annotationfiles to be created e.g. annotation/sample/sample_",
    type=str)
args = parser.parse_args()

dir_out=args.dir_out
if re.search("/",dir_out) and not os.path.isdir(dir_out.rsplit('/', 1)[0]):
    os.makedirs(dir_out.rsplit('/', 1)[0])

dmnd_out = args.file_in
pre = args.path_annot
################################################################################

# paths:
ko_links = pre + "ko.links"
ko_names = pre + "konames.txt"
ko_enz_acc = pre + "ko_enzyme.list"
ko_cazy = pre + "ko_cazy.list"
ko_pathway = pre + "ko_pathway.list"

################################################################################
# Functions:
################################################################################
# Create dictionary so files are read only once and not for every hit
# overwrites previous, so make list; first order on ko-nr!!
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

# if a list exist, append items in list + new item as nested list
# from_list = [a,b,c,a] ; item = A ; new_list = [ [b,E],[b,A] ]
# ===> new_list = [ [b,E],[b,A],[a,A],[c,A] ]
def add_if_new(from_list, list_to_add, new_list):
    if from_list:
        for _ in from_list:
            check=list_to_add.copy() # add copy(), because otherwise original changes with reference
            check.append(_)
            if check not in new_list:
                new_list.append(check)

# write list to file, every item being a line, if nested list, items tab seperated
def list_writer(W_list,file):
    fileObj = open(file, mode='w')
    with fileObj:
        file_writer = csv.writer(fileObj,delimiter='\t', lineterminator = '\n')
        for _ in W_list:
            file_writer.writerow(_)
    fileObj.close()



################################################################################

print("Building dictionaries")
D_links = dict_maker(ko_links)
D_names = dict_maker(ko_names)
D_enz_acc = dict_maker(ko_enz_acc)
D_cazy = dict_maker(ko_cazy)
D_path = dict_maker(ko_pathway)

# empty list for writing annotation
# not dictionaries, because keys need to be unique
W_names = []
W_enz_acc = []
W_cazy = []
W_paths = []


print("Creating files...")
# to summarize all above in 1 file
summary_file = open(dir_out + "summary.tsv", mode='w')
smm_writer = csv.writer(summary_file,delimiter='\t', lineterminator = '\n')

c=1
unique=["",[""]] # to have correct structure for first loop
with summary_file as summary:
    # header row for csv
    smm_writer.writerow(["contig_ID","dmnd_ID", "KO_ID", "name", "EC", "cazy","pathways"])
    for entry in dmnd_out.readlines():
        #if c >10:
        #    break
        # let user know script is still running
        if c%200000 == 0:
            print("Processed {} lines".format(c))
        entry = entry.split('\t')
        contig_id=entry[0]
        dmnd_id=entry[1]
        ko_id=D_links.get(dmnd_id)
        # keep only unique combinations, top hit first, only if ko_id found
        if ko_id and (contig_id not in unique and ko_id not in unique[1]):
        #if [contig_id,ko_id] != unique and ko_id: # this will give doubles if multiple ko_ids
        # can entry have multiple ko_id's? # yes
            for ko in ko_id: 
                # name: need ko_id without ko:
                # ko:K11070 --> [ko, K11070]
                name=D_names.get(ko.split(':')[1]) # using get will not give an error if key does not exist
                add_if_new(name,[contig_id,ko],W_names)

                enz_acc=D_enz_acc.get(ko)
                add_if_new(enz_acc, [contig_id], W_enz_acc)

                cazy=D_cazy.get(ko)
                add_if_new(cazy, [contig_id], W_cazy)

                paths=D_path.get(ko)
                # keep only maps
                maps=paths # for uniformity, to not attach empty list if should be 'None' 
                if paths:
                    maps=[]
                    for _ in paths:
                        if re.search("map",_):
                            maps.append(_)
                add_if_new(maps, [contig_id], W_paths)

                # keeping prefix ko:, path: ... for readibility in large files?
                smm_writer.writerow([contig_id,dmnd_id,ko,name,enz_acc,cazy,maps])
                #post_data = {'contig_id' : contig_id,'dmnd_id' : dmnd_id,'ko_id' : ko,'name' : name,'enz_acc' : enz_acc,'cazy' : cazy,'pathways' : maps}
                #info.insert_one(post_data)
            unique = [contig_id,ko_id]
        c+=1
    dmnd_out.close()
# Close tsv files
summary.close()

# to create files containing 2 columns, 1st= contig_id , 2nd= kegg ref
# check for duplicates!!
list_writer(W_names, dir_out + "names")
list_writer(W_enz_acc, dir_out + "enz_acc")
list_writer(W_cazy, dir_out + "cazy")
list_writer(W_paths, dir_out + "pathways")

print("Processed {} lines in total.".format(c))

