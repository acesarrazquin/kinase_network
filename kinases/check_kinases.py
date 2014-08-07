#!/usr/bin/python

#############################################################################
#
# Check which protein kinases from kinome.com are not included in UniProt's
# dataset, and why (manually check those cases)
#
# Input: uniprot_kinase_list and kinase_list from kinome.com (query all human
#	kinases and get table)
#############################################################################

from drugs import *


name2uniprotid = getName2UniprotDic()
#uniprotid2name = getUniprot2NameDic()

# Read kinome file
count = 0
index = 0
print('\nNot found kinases in UniProt kinase file:')
with open('kinases/kinome_kinases_20140218.txt') as kinfile:
    for line in kinfile:
        index += 1
        split_line = line.split('\r\n')[0].split('\t')
        prot_id = split_line[0]

        synonyms = split_line[2]
        synonyms = re.sub('[ "]', '', synonyms)
        synonyms_list = list(set(filter(None, synonyms.split(',') + [prot_id])))

        flag = 0
        for synonym in synonyms_list:
            if synonym.upper() in name2uniprotid.keys():
                flag = 1
                break

        if flag == 0:
            count += 1
            print('\t' + prot_id + ', synonyms: ' + ','.join(synonyms_list))

print('\nTotal not found: ' + str(count) + ' of ' + str(index) + '\n')
