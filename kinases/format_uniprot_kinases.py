#!usr/bin/python

import sys
import re

file_name = sys.argv[1]

infile = open(file_name, 'r')
all_lines = infile.read()

split_lines = all_lines.split('\r')

kinases_types = {} #dictionary for types and descriptions
kinases_list = []

atypical = 0
#type_pattern = '\*'
group_pattern = '===='
protein_pattern = '_HUMAN'

for i in range(0, len(split_lines)):
    line = split_lines[i]

#    if re.search('Atypical', line):
#        atypical =1
#   
#    if re.search(type_pattern, line):
#        kinase_type_description = line.rstrip().split('* ')[1].split('; ')
#        if atypical == 1:
#            kinase_type = 'Atypical: ' + kinase_type_description[0]
#        else:
#            kinase_type = kinase_type_description[0]
#        if len(kinase_type_description) > 1:
#            kinase_description = kinase_type_description[1]
#        else:
#            kinase_description = ''
#
#        kinases_types[kinase_type] = kinase_description

    if re.search(group_pattern, line) and re.search(group_pattern, split_lines[i+2]):
        kinase_type = re.split('\\n| protein| kinase| family| Ser/Thr', split_lines[i+1])[1]
        print(kinase_type)    
        for j in range(i+3, len(split_lines)):
            line2 = split_lines[j]
            if re.search(protein_pattern, line2):
                split_kinase = line2.split()
                kinase_name = split_kinase[0]
                kinase_uniprotid = split_kinase[2][1:len(split_kinase[2])-1]
                
                kinase = [kinase_name, kinase_uniprotid, kinase_type]
                kinases_list.append(kinase)
            if re.search(group_pattern, line2):
                break
infile.close()

outfile_name = file_name.split('.txt')[1]+'_formated.txt'
outfile = open(outfile_name, 'w')
outfile.write('Name'+'\t'+'UniProt_ID'+'\t'+'Type'+'\n')
for line in kinases_list:
    outfile.write(line[0]+'\t'+line[1]+'\t'+line[2]+'\n')
outfile.close()    
