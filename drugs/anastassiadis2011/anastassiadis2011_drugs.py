#!/usr/bin/python
##########################################################################
#   This file extracts the kinase - drug interactions from the Anastassiadis
#   2011 paper.
#
#   Input: automatic
#	   1. Supplementary table from the paper (tab separated)
#          2. Table of kinases (tab separated), extracted from the PDF of
#             supplementary material
#          3. File of kinases from UniProt, including synonyms etc.
#
#   Output: list of drug to kinase interactions, tab separated
#
#   Some remarks:
#       - names are first mapped to the HUGO nomenclature using the tables
#           from the article, then to UniProt using the kinases file
#       - drug identifiers are mapped to PubChem ids using the CAS numbers
#           as name synonyms (Entrez utilities)
##########################################################################
import sys
import re
from drugs import *
import datetime

def getName2StdNameANDUniprotDic(anastas_kinases_filename):
    result = {}
    with open(anastas_kinases_filename) as anaskin_file:
        for line in anaskin_file:
            split_line = line.rstrip().split('\t')
            prot_name = split_line[0].upper().split('/')[0]
            prot_std_name = split_line[1]
            try:
                uniprot_id = split_line[3]
            except:
                uniprot_id = split_line[2]
            result[prot_name] = [prot_std_name, uniprot_id]
            result[prot_std_name] = [prot_std_name, uniprot_id]
    return result

def writeLine2File(outfile):
    outfile.write(compound_name+'\t'+uniprot_id+'\t'+source_name+'\t'+prot_name+'\t'+
                  pmid+'\t'+dates+'\t'+type_+'\t'+pchem_cid+'\t'+target_form+'\t'+
                  strength+'\t'+strength_type+'\t'+strength_direction+'\t'+strong_target+'\t'+
                  technology+'\t'+cell_type+'\t'+sources+'\n')

def writeHeader2File(outfile):
    outfile.write('source'+'\t'+'target'+'\t'+'source_name'+'\t'+'target_name'+'\t'+'PMIDs'+'\t'+
                  'dates'+'\t'+'type'+'\t'+'source_id'+'\t'+'target_form'+'\t'+
                  'strength'+'\t'+'strength_type'+'\t'+'strength_direction'+'\t'+'strong_target'+'\t'+
                  'technology'+'\t'+'cell_type'+'\t'+'sources'+'\n')

#///////////////////////////////// MAIN ///////////////////////////////////
today = datetime.datetime.now().strftime('%Y%m%d')

anastas_filename = 'anastassiadis2011_table3.txt'
anastas_kinases_filename = 'anastassiadis2011_kinases.txt'

outfilename = 'anastassiadis2011_drugs_formatted-'+today+'.txt'

# some parameters for output file
strength_type = '% remaining kinase activity at [drug] = 500nM' # average IC50 for primary targets is 66nM (therefore much excess in this experiment)
strength_direction = 'low'
cell_type = ''
technology = 'catalytic assay'
pmid = '22037377'
sources = pmid
dates = '2011'
type_ = 'DPI'
target_form = ''

# mapping dictionaries
print('\nGetting UniProt mappings...')
name2uniprot_dic = getName2UniprotDic()
uniprot2name_dic = getUniprot2NameDic()
name2stdname_uniprot_dic = getName2StdNameANDUniprotDic(anastas_kinases_filename) # Patch for names: kinase names to HUGO names and uniprot_ids
name2cid_stdname_dic = getDrugName2CIDandStdName_Dic()

# read table from publication and write output file
with open(anastas_filename) as anastasfile, open(outfilename, 'w') as outfile:
    writeHeader2File(outfile)
    lines = anastasfile.readlines()

    header = lines[0].split('\t')
    compound_names = filter(None, header[1:len(header) - 1])

    casline = lines[1].split('\t')
    cas_nbs = filter(None, casline[1:len(casline) - 1])

    # retrieve PuChemCIDs from 'list_drugs.txt' file or Entrez using CAS numbers or otherwise name
    print('\nRetrieving PubChem Compound IDs for drugs...')
    pchem_cids = []
    for i in range(0, len(compound_names)):
        cas = cas_nbs[i]
        compound_name = ' '.join(compound_names[i].split()) # to substitute sequences of spaces for only one space

        if compound_name == 'Alsterpaullone 2-Cyanoethyl':
            compound_name = '"Alsterpaullone, 2-Cyanoethyl"'

        try: # check if the name is in the list
            pchem_cid = name2cid_stdname_dic[compound_names[i]]['cid']
        except:
            pchem_cid = '; '.join(getCIDFromCAS(cas))

            if pchem_cid == '':
                pchem_cid = '; '.join(getCIDFromName(compound_name))

                if pchem_cid == '':
                    pchem_cid = ''
                    print('\t' + compound_name + '(' + cas + ') : no PubChem CID retrieved!')
        pchem_cids.append(pchem_cid)
    # parse file
    print('\nParsing drug interactions file...')
    for line in lines[2:len(lines)]:

        split_line = line.rstrip().split('\t')
        prot_names = re.sub(' ', '', split_line[0])
        prot_name = name2stdname_uniprot_dic[prot_names.upper().split('/')[0]][0].split('/')[0] # because some kinases have two names separated by a /
        try:
            uniprot_id = name2uniprot_dic[prot_name.upper()]
        except:
            try:
                uniprot_id = name2stdname_uniprot_dic[prot_name.upper()][1]
            except:
                print('\t' + prot_names + '(' + prot_name + ') : no UniProt ID retrieved!')
                uniprot_id = ''
                continue

        if uniprot_id != '':
            prot_name = uniprot2name_dic[uniprot_id]

        strengths = split_line[1:len(split_line)]
        for i in range(0, len(strengths)):
            strength = strengths[i]
            if strength and float(strength) <= 50:
                compound_name = ' '.join(compound_names[i].split())
                pchem_cid = pchem_cids[i]
                source_name = compound_name

                if float(strength) <= 25:
                    strong_target = 'yes'
                else:
                    strong_target = 'no'
                writeLine2File(outfile)


