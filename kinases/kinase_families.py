#!/usr/bin/python

##########################################################################################
#
#	Annotate the Kinases from Uniprot with group, family, subfamyly of Kinome.es file
#
#
#########################################################################################
name2all = {}
fields = ['name', 'std_symbol', 'group', 'fam', 'subfam', 'synonyms', 'entrez_id', 'refs']
for field in fields:
    name2all[field] = []
with open('kinases/Kincat_Hsap.08.02-reduced.txt') as kinomefile:
    lines = kinomefile.readlines()
    for line in lines[1:len(lines)]:
        split_line = line.split('\n')[0].split('\t')

        name = split_line[0]
        group = split_line[2]
        fam = split_line[3]
        subfam = split_line[4]
        entrez_id = split_line[5]
        std_symbol = split_line[6]
        synonyms = split_line[7]
        refs = split_line[8]

        for field in fields:
            name2all[field].append(vars()[field])

fields = ['name', 'approved_symbol', 'synonyms_all', 'synonyms_filtered', 'approved_name', 'type_', 'uniprot_id', 'entrez_id', 'hgnc_id']
with open('kinases/uniprot_kin_formatted_annotated.txt') as uniprotfile, open('uniprot_kin_formatted_annot.txt', 'w') as outfile:
    outfile.write('name\tapproved_symbol\tapproved_name\tuni_type\tgroup\tfam\tsubfam\tsynonyms\tuniprot_id\tentrez_id\thgnc_id\trefs\n')
    lines = uniprotfile.readlines()
    for line in lines[1:len(lines)]:
        index = -1
        split_line = line.split('\n')[0].split('\t')
        for i in range(0, len(fields)):
            vars()[fields[i]] = split_line[i]

        # search for protein in the kinase.com dictionary, by name or entrez id

        if name in name2all['name']:
            index = name2all['name'].index(name)
        elif name in name2all['std_symbol']:
            index = name2all['std_symbol'].index(name)
        elif entrez_id != '' and entrez_id in name2all['entrez_id']:
            index = name2all['entrez_id'].index(entrez_id)
        #print(name + ' found by entrez_id' + str(entrez_id) + ' : ' + name2all['std_symbol'][index])
        else:
            for synonym in synonyms_all.split(', '):
                if synonym in name2all['name']:
                    index = name2all['name'].index(synonym)
                    #print(name + ' found with synonym ' + synonym)
                    break

        if index != -1:
            group = name2all['group'][index]
            fam = name2all['fam'][index]
            subfam = name2all['subfam'][index]

            synonyms_kinome = name2all['synonyms'][index].split('|')
            synonyms_uni = synonyms_all.split(', ')
            stdname = name2all['std_symbol'][index]
            kinoname = name2all['name'][index]

            synonyms_list = list(filter(None, set(synonyms_uni+synonyms_kinome+[stdname]+[kinoname]+[name])))
            if '-' in synonyms_list:
                synonyms_list.remove('-')

            synonyms = '|'.join(synonyms_list)
        else:
            print(name + ' not found in kinase.com file')
            group = ''
            fam = ''
            subfam = ''

            synonyms = '|'.join(synonyms_all.split(', '))

        outfile.write(name+'\t'+approved_symbol+'\t'+approved_name+'\t'+
                      type_+'\t'+group+'\t'+fam+'\t'+subfam+'\t'+synonyms+'\t'+
                      uniprot_id+'\t'+entrez_id+'\t'+hgnc_id+'\t'+refs+'\n')









