#!/usr/bin/python

#####################################################################################################
# Script to separate strength values from drug protein interaction files between primary drug targets
# as defined by davis2011 and rest of targets.
#
# Input: - drug protein interactions file (as in "formatted_updated")
#
# Output: - drug protein interactions as: drug \t protein_target \t strength \t ["primary"|"not"]
#####################################################################################################
from drugs import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st

# choose file:
infilename = sys.argv[1]
outfilename = infilename.split('formatted_updated.txt')[0] +'distributions.txt'

# get dictionaries:
name2uniprot = getName2UniprotDic()

## get distributions for primary targets and all targets to choose threshold, based on davis file
##finish this
primary_targets = {}
with open('davis2011_compounds.txt') as compfile:
    lines = compfile.readlines()
    for line in lines[1:len(lines)]:
        split_line = line.split('\n')[0].split('\t')
        if split_line[1] != '':
            compound_name = split_line[1]
        else:
            compound_name = split_line[0]

        if compound_name == 'INCB018424':
            compound_name = 'INCB18424'

        primary_targets[compound_name] = []
        targets = filter(None, [split_line[2]] + [split_line[3]] + [split_line[4]])
        for target in targets: # put original name, and split names that have bars
            if target == 'p38-alpha':
                target = 'p38'

            if re.search('/', target):
                split_targets = target.split('/')

                for i in range(0, len(split_targets)):
                    if i>0:
                        split_target = split_targets[0][0:len(split_targets[0])-1] + split_targets[i]
                    else:
                        split_target = split_targets[0]

                    try:
                        uniprot_id = name2uniprot[split_target.upper()]
                        primary_targets[compound_name].append(uniprot_id)
                    except:
                        print('No uniprot ID found for primary target ' + split_target)
            else:
                try:
                    uniprot_id = name2uniprot[target.upper()]
                    primary_targets[compound_name].append(uniprot_id)
                except:
                    print('No uniprot ID found for primary target ' + target)
                #if target == 'ABL1(T315I)': # I will skip it.
                #	uniprot_id = 'ABL1(T315I)'
                #	print('\t but included')
                #	primary_targets[compound_name].append(uniprot_id)


with open(infilename) as interactions_file, open(outfilename, 'w') as outfile:
    lines = interactions_file.readlines()
    for line in lines[1:len(lines)]:
        split_line = line.split('\n')[0].split('\t')
        kd = split_line[9]
        drug = split_line[0]
        prot = split_line[1]
        prot_name = split_line[3]

        if re.search('davis', infilename):
            form = split_line[8]
        else:
            form = ''

        if form == '':
            try:
                if prot in primary_targets[drug]:
                    outfile.write(drug + '\t' + prot_name + '\t' + str(float(kd)) + '\tprimary\n')
                else:
                    outfile.write(drug + '\t' + prot_name + '\t' + str(float(kd)) + '\tnot\n')
            except:
                outfile.write(drug + '\t' + prot_name + '\t' + str(float(kd)) + '\tnot\n')
            #print('Drug not in primary targets ' + drug)
        else:
            outfile.write(drug + '\t' + prot_name + '\t' + str(float(kd)) + '\tnot\n')

#prim_density = st.gaussian_kde(prim_targets)
#xs = np.linspace(0,1000,10)
#prim_density.covariance_factor = lambda : 0.25
#prim_density._compute_covariance()
#plt.plot(xs,prim_density(xs))
#plt.show()
