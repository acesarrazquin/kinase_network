import re
import sys
import copy


kinases_file_name = sys.argv[1]
genenames_file_name = sys.argv[2]
entrez_file_name = sys.argv[3]
englishwords_file_name = sys.argv[4]


uniprot_attributes_dic = {}
entrez_attributes_dic = {}
fields_list = ['hgnc_id', 'approved_symbol', 'approved_name', 'previous_symbols', 'previous_names', 'synonyms', 'name_synonyms', 'entrez_id', 'uniprot_id']
all_names = [] # list of all names together, to see if they are repeated ones
with open(genenames_file_name) as genenames_file:
    print('\nMapping UniProt identifiers and attributes from HGNC...') 
    for line in genenames_file:
        split_lines = line.split('\r')[0].split('\t')# some lines with empty fields at the end have one less \t. Instead of rstrip(), use split('\r') -- carriage return
        
        if split_lines[0] == 'HGNC ID': # to skip the first line
            continue
        
        for i in range(0,len(fields_list)):
            vars()[fields_list[i]] = split_lines[i]
        
        all_synonyms = filter(None, approved_symbol.split(', ') + previous_symbols.split(', ') + synonyms.split(', ')) #list of synonms
        all_names.extend(all_synonyms)
        uniprot_attributes_dic[uniprot_id] = [hgnc_id, entrez_id, approved_symbol, approved_name, all_synonyms]
        entrez_attributes_dic[entrez_id] = [hgnc_id, entrez_id, approved_symbol, approved_name, all_synonyms] 
del uniprot_attributes_dic[''] # erase entries that are empty or just contain a '-' symbol
del uniprot_attributes_dic['-']    
del entrez_attributes_dic['']


englishwords = []
print('\nPreparing filters...')
with open(englishwords_file_name) as englishwords_file:
    for line in englishwords_file:
        word = line.rstrip().upper()
        englishwords.append(word)

unique_names = list(set(all_names))
repeated_names = []
for name in unique_names:
    if all_names.count(name) > 1:
        repeated_names.append(name)


uniprot_entrez_dic = {}
with open(entrez_file_name) as entrez_file:
    print('\nMapping Entrez Gene (Gene_ID) and UniProt identifiers...')
    for line in entrez_file:
        split_lines = line.rstrip().split('\t')
        uniprot_id = split_lines[0]
        entrez_id = split_lines[2].split('; ') # if more than one entrez id, take all in list. Anyway it is always a list

        uniprot_entrez_dic[uniprot_id] = entrez_id
   


outfile_name = kinases_file_name.split('.txt')[0]+'_annotated2.txt'
all_kinase_names = []
with open(kinases_file_name) as kinases_file, open(outfile_name, 'w') as outfile:

    outfile.write('Name'+'\t'+'ApprovedSymbol'+'\t'+'Synonyms_all'+'\t'+'Synonyms_filtered'+'\t'+'ApprovedName'+'\t'+
                  'Type'+'\t'+'UniProt_ID'+'\t'+'Entrez_ID'+'\t'+'HGNC_ID'+'\n')
    print("\nLooking for UniProt kinases synonyms and filtering...")
    
    for line in kinases_file:
        split_lines = line.rstrip().split('\t')
        if split_lines[0] == 'Name':
            continue # avoid first line
                        
        kinase_name = split_lines[0]
        uniprot_id = split_lines[1]
        kinase_type = split_lines[2]
        entrez_ids = uniprot_entrez_dic[uniprot_id] 

        try:
            attributes = uniprot_attributes_dic[uniprot_id]
        except:
            print('\t'+kinase_name + ': not found!')
            for entrez_id in entrez_ids:
                try:
                    attributes = entrez_attributes_dic[entrez_id]
                    print('\t\tbut recovered from Entrez_ID!')
                    break
                except:
                    attributes = ['' ,'', '', '', []]
        
        hgnc_id = attributes[0]
        entrez_id = attributes[1]
        if entrez_id == '':
            entrez_id = ', '.join(entrez_ids)
        else:
            if len(entrez_ids) > 2 and entrez_id in entrez_ids:
                entrez_id = ', '.join(entrez_ids)
        approved_symbol = attributes[2]
        approved_name = attributes[3]
        all_synonyms = attributes[4]
        
        # apply filters:
        filt_synonyms = copy.copy(all_synonyms)
        for synonym in filt_synonyms:
            synonym_ = synonym.upper()
            if synonym in repeated_names:
                print(kinase_name + ' : ' + synonym + ' removed -- repeated')
                filt_synonyms.remove(synonym)
            elif synonym_ in englishwords:
                print(kinase_name + ' : ' + synonym + ' removed -- English word')
                filt_synonyms.remove(synonym)

        # add the initial name if not present        
        if kinase_name not in all_synonyms:
            all_synonyms.insert(0, kinase_name)
        if kinase_name not in filt_synonyms:
            filt_synonyms.insert(0, kinase_name)
        allsynonyms = ', '.join(all_synonyms)
        filtsynonyms = ', '.join(filt_synonyms)
        all_kinase_names.append(kinase_name)

        outfile.write(str(kinase_name)+'\t'+str(approved_symbol)+'\t'+str(allsynonyms)+'\t'+str(filtsynonyms)+'\t'+
                      str(approved_name)+'\t'+str(kinase_type)+'\t'+str(uniprot_id)+'\t'+str(entrez_id)+'\t'+str(hgnc_id)+'\n')

# How many of the kinase names would have been removed?
print('\nMain kinase names that would have been removed if filteres were applied:')
for name in all_kinase_names:     
    if name in repeated_names and name in englishwords:
        print('\t'+name+': repeated and English word!')
    elif name in repeated_names:
        print('\t'+name+': repeated!')
    elif name in englishwords:
        print('\t'+name+': English word!')                                                                                                                                                   
