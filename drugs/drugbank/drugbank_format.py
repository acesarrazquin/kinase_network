#!/usr/bin/python
import sys
import re

selected_features = ['DrugBank_ID', 'Generic_Name', 'Synonyms', 'Brand_Names', 'Chemical_IUPAC_Name', 'Chemical_Formula', 'Drug_Category', 
                     'Drug_Type', 'Protein_Binding', 'Molecular_Weight_Avg', 'Drug_Reference', 
                     'Primary_Accession_No', 'Secondary_Accession_No', 'ChEBI_ID', 'CAS_Registry_Number', 'PubChem_Compound_ID', 
                     'PubChem_Substance_ID', 'InChI_Identifier', 'InChI_Key', 'PharmGKB_ID', 'LIMS_Drug_ID', 'DPD_Drug_ID_Number', 
                     'GenBank_ID', 'HET_ID', 'KEGG_Compound_ID', 'KEGG_Drug_ID', 'PDB_Experimental_ID', 'PDB_Homology_ID', 
                     'SwissProt_ID', 'SwissProt_Name', 'Description', 'Mechanism_Of_Action'] 

outfilename = 'drugbank_formatted.txt'
with open('drugbank.txt') as drugbankfile, open(outfilename, 'w') as outfile:
    
    outfile.write('\t'.join(selected_features) + '\n')
    
    
    drugbank = drugbankfile.read()
    split_file = drugbank.split('#BEGIN_DRUGCARD ')


    del split_file[0]
    
    for drugcard in split_file:

        split_features = drugcard.split('#')
        db_id = split_features.pop(0).split('\n')[0]
        del split_features[len(split_features)-1]

        print(db_id)
        drug_features_dic = {}
        drug_features_dic['DrugBank_ID'] = db_id
        for feature in split_features:
            split_lines = filter(None, feature.split('\n'))
             
            feature_name = split_lines.pop(0)[1:-1]
            
            if feature_name == 'Drug_Reference':
                pmids = []
                for split_line in split_lines:
                    split_fields = split_line.split('\t')
                    if re.match('[0-9]+', split_fields[0]):
                        pmids.append(split_fields[0])
                value = '; '.join(pmids)
                drug_features_dic[feature_name] = value

            elif feature_name in selected_features:
                value = re.sub('\t', '', '; '.join(split_lines))
                
                drug_features_dic[feature_name] = value 
        
        for selected_feature in selected_features:
            outfile.write(drug_features_dic[selected_feature] + '\t')
        outfile.write('\n')

