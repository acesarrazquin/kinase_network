#!/usr/bin/python

#####################################################################################################
# File to extract drug2protein interactions from the bantscheff2007 paper.
#	- only for bosutinib, dasatinib and imatinib (the other experiments, that included more drugs,
#	  are less reliable...) -> PubChem CIDs added manually
#	- filter out the IC50s below 1uM
#	- map kinase names to UniProt IDs using the initial kinase list from UniProt
#	
# Input: automatic (two supplementary tables extracted from the pdf of bantscheff2007 -lysate and cell-
#	 culture - and kinases file)
#
# Output: file of interactions with one protein-drug interaction per row, tab separated
#
#####################################################################################################
import datetime
import sys
import re
from drugs import *

def writeLine2File(outfile):
    outfile.write(compound_name+'\t'+uniprot_id+'\t'+source_name+'\t'+kin_name+'\t'+
                  pmid+'\t'+dates+'\t'+type_+'\t'+pchem_cid+'\t'+target_form+'\t'+strength+'\t'+strength_type+'\t'+
                  strength_direction+'\t'+strong_target+'\t'+technology+'\t'+cell_type+'\t'+sources+'\n')

def writeHeader2File(outfile):
   outfile.write('source'+'\t'+'target'+'\t'+'source_name'+'\t'+'target_name'+'\t'+'PMIDs'+'\t'+
                 'dates'+'\t'+'type'+'\t'+'source_id'+'\t'+'target_form'+'\t'+
				 'strength'+'\t'+'strength_type'+'\t'+'strength_direction'+'\t'+'strong_target'+'\t'+
				 'technology'+'\t'+'cell_type'+'\t'+'sources'+'\n')

#///////////////////////////// MAIN //////////////////////////////////////////////////////////////////
today = datetime.datetime.now().strftime('%Y%m%d')

lysate_filename = 'bantscheff2007-tableS4.txt'
culture_filename = 'bantscheff2007-tableS5.txt'
outfilename = 'bantscheff2007_drugs_formatted-'+today+'.txt'

# constant fields for output file:
strength_type = 'IC50[nM]' # concentration of drug at which half maximal competition of kinobead is observed
strength_direction = 'low'
pmid = '17721511'
sources = pmid
dates = '2007'
type_ = 'DPI' # ?
cell_type = 'K562'


print('\nGetting UniProt ID mappings...')
# mapping dictionaries
name2uniprot_dic = getName2UniprotDic()
uniprot2name_dic = getUniprot2NameDic()

drugs = ['Bosutinib', 'Dasatinib', 'Imatinib']
pubchem_cids = ['5328940', '3062316', '5291']
with open(outfilename, 'w') as outfile:
    writeHeader2File(outfile)
    for infilename in (lysate_filename, culture_filename):

        with open(infilename) as infile:
            if infilename == lysate_filename:
                technology = 'chemprot - kinobeads competition assay - lysate'
            else:
                technology = 'chemprot - kinobeads competition assay - cell culture'
            
            print('\nParsing drug interaction file ' + infilename + '...') 		
	    for line in infile:
                split_line = line.rstrip().split('\t')
               
                try:
                    kinase_type = split_line[2]
                    prot_name = split_line[1].upper()
                except:
                    continue

                if not re.match('IPI', split_line[0]) and prot_name != 'BCR-ABL':
                    continue
                if re.match('SIMILAR', split_line[1]):
                    continue

                if kinase_type or prot_name == 'BCR-ABL': # here I take BCR/ABL...as it's not listed as kinase in the file. But I wont give accesion nb
                    alternative_name = ''
                    target_form = ''
                    Bosutinib_score = split_line[14]
                    Dasatinib_score = split_line[25]
                    Imatinib_score = split_line[34]
            
                    try:
                        uniprot_id = name2uniprot_dic[prot_name]
                    except:
                        if prot_name == 'ZAK':
                            uniprot_id = 'Q9NYL2'
                            print('\t' + prot_name + ': UniProt ID not found - but recovered: synonym of MLTK')
                        elif prot_name == 'BCR-ABL':
                            uniprot_id = 'P00519'
                            target_form = 'BCR-ABL'
                            print('\t' + prot_name + ': included using ABL1 UniProt ID')
                        else:
                            print('\t' + prot_name + ': UniProt ID not found')
                            continue # i dont include the kinases not in the original list from uniprot
		    		
                    try:
                        kin_name = uniprot2name_dic[uniprot_id]
                    except:
                        kin_name = prot_name

                    for i in range(0, len(drugs)):
                        score = vars()[drugs[i]+'_score']
                        if score != '>10' and score != '>5' and float(score) <= 1: # only consider IC50 <= 1uM
                            compound_name = drugs[i]
                            source_name = compound_name
                            pchem_cid = 'CID:' + pubchem_cids[i]
                            
                            if float(score) <= 0.2:
                                strong_target = 'yes'
                            else:
                                strong_target = 'no'
                            strength = str(int(float(score)*1000)) # convert uM to nM
                            

                            writeLine2File(outfile)

# Maybe I could parse this other file with more drugs (but the data are less reliable)
