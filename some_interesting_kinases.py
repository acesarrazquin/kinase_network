#!/usr/bin/python
import sys
from drugs import *

mode = sys.argv[1]
# dictionary and list of kinases
id2name = getUniprot2NameDic()


# dictionary of interactions
fields = ['source', 'target', 'source_name', 'target_name', 'PMIDs', 'dates', 'type_']
interactions = {}
for field in fields:
    interactions[field] = []
with open('full-network-20140130.txt') as infile:
    lines = infile.readlines()

    for line in lines[1:len(lines)]:
        split_line = line.split('\n')[0].split('\t')

        for i in range(0, len(fields)):
            interactions[fields[i]].append(split_line[i])

kinases = id2name.keys()
if mode == 'PPI':
    # 10 kinases most connected by PPI to other kinases:

    kin_counts = {}
    kin_interactors = {}
    for kinase in kinases:
        count = 0
        interactors = []
        for i in range(0, len(interactions['source'])):
            if interactions['type_'][i] == 'PPI':# change this if necessary
                if interactions['source'][i] == kinase and interactions['target'][i] in kinases:
                    count += 1
                    interactors.append(interactions['target'][i])
                elif interactions['target'][i] == kinase and interactions['source'][i] in kinases: #for PPI, not for KSI
                    count += 1
                    interactors.append(interactions['source'][i])

        kin_counts[kinase] = count
        kin_interactors[kinase] = interactors

    best = sorted(kin_counts.items(), key = lambda x : x[1], reverse=True)[0:10]
    print('\nKinases with most PPI to other kinases:')
    for i in best:
        name = id2name[i[0]]
        print(name + ' (' + i[0] + ') : ' + str(i[1]))
        for j in kin_interactors[i[0]]:
            name2 = id2name[j]
            print('\t' + name2 + ' (' + j+ ')')

elif mode == 'nonkin':
    nonkin_ksi = {}
    for i in range(0, len(interactions['source'])):
        if interactions['type_'][i] == 'PPI':# or interactions['type_'][i] == 'computeKSI':
            if interactions['source'][i] in kinases and interactions['target'][i] not in kinases:
                nonkin = interactions['target'][i]
                kin = interactions['source'][i]
                if nonkin in nonkin_ksi.keys():
                    if kin not in nonkin_ksi[nonkin]:
                        nonkin_ksi[nonkin].append(kin)
                else:
                    nonkin_ksi[nonkin] = [kin]
            elif interactions['source'][i] not in kinases and interactions['target'][i] in kinases:
                nonkin = interactions['source'][i]
                kin = interactions['target'][i]
                if nonkin in nonkin_ksi.keys():
                    if kin not in nonkin_ksi[nonkin]:
                        nonkin_ksi[nonkin].append(kin)
                else:
                    nonkin_ksi[nonkin] = [kin]
    best = sorted(nonkin_ksi.items(), key = lambda x : len(x[1]), reverse=True)[0:10]
    print('\nNon-kinases most targetted by kinases (PPI):')
    for i in best:
        name = i[0]
        print('\t' + name + ': ' + str(len(i[1])))
    #for j in nonkin_ksi[name]:
    #	name2 = id2name[j]
    #	print('\t' + name2 + ' (' + j+ ')')
elif mode == 'drugs':
    kindrugs = {}
    for i in range(0, len(interactions['source'])):
        if interactions['type_'][i] == 'DPI':
            drug = interactions['source'][i]
            kin = interactions['target'][i]
            if kin in kindrugs.keys():
                kindrugs[kin].append(drug)
            else:
                kindrugs[kin] = [drug]

    best = sorted(kindrugs.items(), key = lambda x : len(x[1]), reverse=True)[0:10]
    print('\nKinases most targetted by drugs (DPI):')
    for i in best:
        name = id2name[i[0]]
        print(name + ' (' + i[0] + ') : ' + str(len(i[1])))
        for j in kindrugs[i[0]]:
            print('\t' + j)

