#!/usr/bin/python
###########################################################################################################
#
#
#
#
#
#
###########################################################################################################
import re
import sys

def writeLine2File(outfile):
    outfile.write(compound_name+'\t'+kinase_id+'\t'+source_name+'\t'+kinase_name+'\t'+
	              pmid+'\t'+dates+'\t'+type_+'\t'+pchem_cid+'\t'+target_form+'\t'+
				  strength+'\t'+strength_type+'\t'+strength_direction+'\t'+technology+'\t'+cell_type+'\t'+sources+'\n')

def writeHeader2File(outfile):
	outfile.write('source'+'\t'+'target'+'\t'+'source_name'+'\t'+'target_name'+'\t'+'PMIDs'+'\t'+
				  'dates'+'\t'+'type'+'\t'+'source_id'+'\t'+'target_form'+'\t'+
				  'strength'+'\t'+'strength_type'+'\t'+'strength_direction'+'\t'+'technology'+'\t'+'cell_type'+'\t'+'sources'+'\n')

#//////////////////////////////////////// MAIN //////////////////////////////////////////////////////////////
import sys
import re
from drugs import *
import datetime

today = datetime.datetime.now().strftime('%Y%m%d')

targets_filename = 'drugbank/all_target_ids_with_known_action.csv' # -- list all targets different IDs, and the drugs associated with them
outfilename = 'drugbank_kinasedrugs_formatted_filt-'+today+'.txt'


# parameters:
sources = 'drugbank'
cell_type = ''
technology = ''
strength_direction = ''
strength_type = ''
strength = ''
target_form = ''
type_ = 'DPI'
dates = ''
pmid = ''


# dictionaries:
drugBank_dic = getDrugBankDic()
name2uniprot_dic = getName2UniprotDic()
uniprot2name_dic = getUniprot2NameDic()

# reading targets from drugs in DrugBank
fields = ['ID', 'Name', 'Gene_Name', 'GenBank_Protein_ID', 'GenBank_Gene_ID', 'UniProt_ID', 'Uniprot_Title',
          'PDB_ID', 'GeneCard_ID', 'GenAtlas_ID', 'HGNC_ID', 'HPRD_ID', 'Species_Category', 'Species', 'Drug_IDs']
uniprot_drug_dic = {} #maps uniprot identifier to the list of drugs that target it
print('\nLooking at targets from drugs in DrugBank...')
with open(targets_filename) as targets_file:
    for line in targets_file:
        if re.match('ID', line):
            continue
        split_lines = line.split('\n')[0].split(',') # some fields contain names with commas! but they are quoted
        new_fields = []
        for split_line in split_lines:
            if re.search('^\"', split_line):
                first = split_line[1:len(split_line)]
                if re.search('\"$', split_line):
                    new_field = first[0:len(first)-1]
                    new_fields.append(new_field)        
            elif re.search('\"$', split_line):
                last = split_line[0:len(split_line)-1]
                new_field = first+last
                new_fields.append(new_field)
            else:
                new_field = split_line
                new_fields.append(new_field)
                                 
        for i in range(0, len(fields)):
            vars()[fields[i]] = new_fields[i]

        drugs = Drug_IDs.split('; ') # list of drugsdd
        uniprot_drug_dic[UniProt_ID] = drugs

if '' in uniprot_drug_dic.keys():
    del uniprot_drug_dic['']



found = 0
all_kinase_drugs = []
drug_kinases_dic = {} # dictionary mapping each drug to the kinases targetted by it
print('\nMapping protein kinases targeted by drugs in DrugBank...')
for uniprot_id in uniprot2name_dic.keys():     
	name = uniprot2name_dic[uniprot_id]
	try:	
		kinase_drugs = uniprot_drug_dic[uniprot_id]
		print('\t' + name + ' (' + uniprot_id + ') found in DrugBank targets')
		found += 1
	except:
		#print('\t' + Name + ' (' + UniProt_ID + ') not in DrugBank targets')
		continue
    
	all_kinase_drugs.extend(kinase_drugs)    
        
	for kinase_drug in kinase_drugs:
		if kinase_drug in drug_kinases_dic.keys():
			drug_kinases_dic[kinase_drug].append([uniprot_id, name])
		else:
			drug_kinases_dic[kinase_drug] = []
			drug_kinases_dic[kinase_drug].append([uniprot_id, name])    

print('\n\tTotal numbef of kinases in DB targets: ' + str(found))
unique_kinase_drugs = list(set(all_kinase_drugs))



all_targetted_kinases = []
print('\nRetrieving and writing protein kinase-targetting drug data...')
with open(outfilename, 'w') as outfile:
	writeHeader2File(outfile)
	length = len(drugBank_dic['Primary_Accession_No'])
	for i in range(0, length):
		#fields = ['ChEBI_ID', 'PubChem_Compound_ID', 'Generic_Name', 'Primary_Accession_No', 'Brand_Names', 'Synonyms',] # many other fields possible
		drug_id = drugBank_dic['Primary_Accession_No'][i] 
		if drug_id not in unique_kinase_drugs: 
			continue
		elif not re.search('kinase inhibitor', drugBank_dic['Description'][i]) and \
			 not re.search('kinase inhibitor', drugBank_dic['Mechanism_Of_Action'][i]):
			continue
		else:	
			kinase_ids_names = drug_kinases_dic[drug_id]
			#chebi_id = drugBank_dic['ChEBI_ID'][i]
			compound_name = drugBank_dic['Generic_Name'][i]
			#if re.search(',', compound_name) and not re.match('"', compound_name):
			#	compound_name = '"' + compound_name + '"'
			compound_name = re.sub(',', '', compound_name)
			source_name = compound_name
			
			pchem_cid = drugBank_dic['PubChem_Compound_ID'][i]
			if pchem_cid == 'Not Available':
				pchem_cids = getCIDFromName(compound_name)
				pchem_cid = ', '.join(pchem_cids)
				if pchem_cid == '':
					print('\n' + compound_name + ': no pubchem CID found!')
					chebi_id = drugBank_dic['ChEBI_ID'][i]
					pchem_cid = drug_id# if no pubmed id, then take any other id

			for kinase_id_name in kinase_ids_names:
				kinase_id = kinase_id_name[0]
				kinase_name = uniprot2name_dic[kinase_id]
			
				all_targetted_kinases.append(kinase_name)
         
				writeLine2File(outfile)

unique_targetted_kinases = list(set(all_targetted_kinases))                                  
print('\nSUMMARY:\n\t'+str(len(unique_kinase_drugs))+' drugs, targetting '+str(len(unique_targetted_kinases))+' protein kinases')
