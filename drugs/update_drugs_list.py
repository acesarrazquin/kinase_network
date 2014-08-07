#!/cm/shared/apps/python/2.7.6/bin/python
####################################################################################
# Update/add new drugs to the drug_list file, which contains a list of drugs with
# several identifiers.
#
# Input: a file with a new list of drugs where: Name_drug \t PChemCID(with header)	
#
# Output: append new data to drug_list file.Fields:
#	  compound_name -> first element of SynonymList from PubChem Compound
#	  pubchem_cid
#	  chembl_id, chebi_id, drugbank_id, inchi, inchikey, iupac_name -> retrieve 
#		from SynonymList using BioPython.Entrez.efetch
#	  synonyms_list -> include the names in the original article if they are
#			   different than name_drug
#
# Remarks:
#  - If some of these identifiers are more than one, separate them with commas
#	
####################################################################################

import datetime
from drugs import *


infilename = sys.argv[1]
listfilename = sys.argv[2]
param = sys.argv[3] # 'article', 'pubchem' -> defines what to take as main name for the drug

today = datetime.datetime.now().strftime('%Y%m%d')
outfilename = 'list_drugs-'+today+'.txt'

fields = ['compound_name', 'pubchem_cid', 'synonym_list', 'chembl_id', 'chebi_id', 'drugbank_id', 'inchi', 'inchikey', 'iupac_name']
with open(infilename) as infile, open(listfilename) as listfile:
    print('\nReading original drug_list file...')
    cid2vals_dic = {}
    lines = listfile.readlines()
    for line in lines[1:len(lines)]: # get the drugs and ids in the current file HEADER
        split_line = line.split('\n')[0].split('\t')
        val_list = []
        for i in range(0, len(fields)):
            vars()[fields[i]] = split_line[i]
            val_list.append(vars()[fields[i]])
        cid2vals_dic[pubchem_cid] = val_list

    print('\nAdding new drugs...')
    lines = infile.readlines()
    for line in lines: # get the new drugs NO HEADER
        split_line = line.split('\n')[0].split('\t')

        orig_name = split_line[0]
        pubchem_cid = split_line[1]
        if re.search(',', orig_name) and orig_name[0] != '"': # in case there are commas in the original name
            orig_name = '"' + orig_name + '"'

        if pubchem_cid in cid2vals_dic.keys(): #if drug already in drug_list, only append name to synonyms if it's different than the compound_name
            database_name = cid2vals_dic[pubchem_cid][0]

            if orig_name.upper() != database_name.upper():
                database_synonyms = cid2vals_dic[pubchem_cid][2]
                if database_synonyms != '' and not re.search(orig_name, database_synonyms):
                    cid2vals_dic[pubchem_cid][2] += '; ' + orig_name
                else:
                    cid2vals_dic[pubchem_cid][2] += orig_name
            print('\tAlready in DB: ' + orig_name +' (' + str(pubchem_cid) +')')
        else: # call Entrez to gather information
            try:
                info = getInfoFromCid(pubchem_cid)
                pchem_name = info['SynonymList'][0]
            except:
                info = {'SynonymList':[], 'InChI':'', 'InChIKey' : '', 'IUPACName':''} #in case no pchemcid...cause the query will form but the substring will fail
                pchem_name = ''
            chembl_ids = []
            chebi_ids = []
            drugbank_ids = []
            for synonym in info['SynonymList']:
                if re.match('CHEMBL[0-9]+', synonym):
                    chembl_ids.append(synonym)
                elif re.match('CHEBI:[0-9]+', synonym):
                    chebi_ids.append(synonym)
                elif re.match('DB[0-9]{5}', synonym):
                    drugbank_ids.append(synonym)
            chembl_id = '; '.join(chembl_ids)
            chebi_id = '; '.join(chebi_ids)
            drugbank_id = '; '.join(drugbank_ids)
            inchi = info['InChI']
            inchikey = info['InChIKey']
            iupac_name = info['IUPACName']

            #pchem_name = info['SynonymList'][0]
            if re.search(',', pchem_name):# in case there is a comma, apply quotes to the whole string
                pchem_name = '"' + pchem_name + '"'
            #don't include the synonym if it's one of those identifiers
            if re.search('CHEMBL[0-9]+', pchem_name) \
                    or re.search('CHEBI:[0-9]+', pchem_name) \
                    or re.search('DB[0-9]{5}', pchem_name):
                pchem_name = ''

            if orig_name.upper() != pchem_name.upper():
                if param == 'pubchem':
                    compound_name = pchem_name
                    synonym_list = orig_name
                elif param == 'article':
                    compound_name = orig_name
                    synonym_list =	pchem_name
            else:
                compound_name = orig_name
                synonym_list = ''

            print('\tAdded: ' + compound_name +' (' + str(pubchem_cid) +')')
            cid2vals_dic[pubchem_cid] = [compound_name, pubchem_cid, synonym_list, chembl_id, chebi_id, drugbank_id,
                                         inchi, inchikey, iupac_name]

with open(outfilename, 'w') as outfile:
    print('\nUpdating drug_list file...')
    for field in fields: # write the header
        outfile.write(field + '\t')
    outfile.write('\n')
    for key in cid2vals_dic.keys(): # write all lines back to drug_list file
        for val in cid2vals_dic[key]:
            outfile.write(str(val) + '\t')
        outfile.write('\n')

