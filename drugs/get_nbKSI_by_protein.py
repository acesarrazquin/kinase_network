#!/usr/bin/python
import sys
from drugs import *


# dictionary of interactions
fields = ['source', 'target', 'source_name', 'target_name', 'PMIDs', 'dates', 'type_']
target2kinase = {}
target2count = {}
id2name = {}
with open('full-network-20140211.txt') as infile:
    lines = infile.readlines()

    for line in lines[1:len(lines)]:
        split_line = line.split('\n')[0].split('\t')

        for i in range(0, len(fields)):
            vars()[fields[i]] = split_line[i]

        if type_ == 'KSI':# or type_ == 'computeKSI':
            if target not in target2kinase.keys():
                target2kinase[target] = [source]
                target2count[target] = 1

            else:
                target2kinase[target].append(source)
                target2count[target] += 1

            if source not in id2name.keys():
                id2name[source] = source_name
            if target not in id2name.keys():
                id2name[target] = target_name

sorted_target2count = sorted(target2count.items(), key = lambda x : x[1], reverse=True)

with open('nb_target2kinase_onlyKSI.txt', 'w') as outfile:
    for i in sorted_target2count:
        outfile.write(i[0]+'\t'+id2name[i[0]]+'\t'+str(i[1])+'\n')
