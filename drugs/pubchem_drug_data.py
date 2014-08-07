#!/usr/bin/python
##################################################################
#
# Script to retrieve data compound CIDs (PubChemCompound IDs)
# from PubChem using the Biopython.Entrez module.
#
# Input:
#   - compounds in a tab separated file: compound \t CID
#
# Output:
#   - completed table with e.g. SMILEs, for further fingerprint
#     extraction and clustering
#
##################################################################
import sys
from drugs import *

infilename = sys.argv[1]
outfilename = infilename.split('.')[0] + '_completed.txt'

fields = ['IUPACName', 'SynonymList', 'CHEMBL_ID', 'CHEBI_ID', 'InChI', 'InChIKey', 'MolecularWeight', 
          'Complexity', 'TotalFormalCharge', 'MolecularFormula', 'CanonicalSmiles', 'IsomericSmiles'] 
           
with open(infilename) as infile, open(outfilename, 'w') as outfile:
    inlines = infile.readlines()
    
    outfile.write('CompoundName\tPubChemCID\t') 
    for field_name in fields:
        outfile.write(field_name + '\t')
    outfile.write('\n')    
    
    for line in inlines[1:len(inlines)]:
        split_line = line.rstrip().split('\t')
        print(split_line)
        compound_name = split_line[0]
        if re.search('derivative', compound_name):
            continue
        pchem_cid = str(split_line[1])
       
        cid_info = getInfoFromCid(pchem_cid)
        
        outfile.write(compound_name + '\t' + str(pchem_cid) + '\t')
        for field in fields:
            if field == 'CHEMBL_ID' or field == 'CHEBI_ID':
                value = []
                for synonym in cid_info['SynonymList']:
                    if re.match('CHEMBL[0-9]+', synonym):
                        value = synonym
                    elif re.match('CHEBI:[0-9]+]', synonym):
                        value = synonym                           
            else:
                value = cid_info[field]
            
            if type(value) == list:
                value = '; '.join(value)
            outfile.write(str(value) + '\t')
        outfile.write('\n')
        
                                                                                                                                     
        
