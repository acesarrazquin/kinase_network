#!/cm/shared/apps/python/2.7.6/bin/python

#########################################################################################
#
# Script to change the names and cids of drug-protein interactions extracted from a file
# according to what is defined in the list_drugs.txt file (repository of all drugs seen
# so far)
#
# Input:- file to modify with standard fields (source, target, source.name, target.name...)	
#		- list_drugs.txt file (automatically opened)
#
# Output: same file but updated
#
#########################################################################################

def writeLine2File(outfile): ##CHECK ORDER OF FIELDS
    outfile.write(compound_name+'\t'+target+'\t'+source_name+'\t'+target_name+'\t'+
                  pmids+'\t'+dates+'\t'+sources+'\t'+type_+'\t'+
                  source_ids+'\t'+target_form+'\t'+strength+'\t'+strength_type+'\t'+
                  strength_direction+'\t'+strong_target+'\t'+technology+'\t'+cell_type+'\n')

def writeHeader2File(outfile):
    outfile.write('source'+'\t'+'target'+'\t'+'source_name'+'\t'+'target_name'+'\t'+'PMIDs'+'\t'+
                  'dates'+'\t'+'sources'+'\t'+'type'+'\t'+
                  'source_IDs'+'\t'+'target_form'+'\t'+'strength'+'\t'+'strength_type'+'\t'+
                  'strength_direction'+'\t'+'strong_target'+'\t'+'technology'+'\t'+'cell_type'+'\n')

#////////////////////////////////////// MAIN /////////////////////////////////////////////

from drugs import *

infilename = sys.argv[1]
drugs_filename = sys.argv[2]
outfilename = infilename.split('.')[0] + '_updated.txt'

fields = ['source', 'target', 'source_name', 'target_name', 'pmids',
          'dates',  'sources', 'type_',
          'source_id', 'target_form', 'strength', 'strength_type',
          'strength_direction', 'strong_target', 'technology','cell_type']

# DrugBank dictionary:
drugbank_dic = getDrugBankDic()
drugname2dbid_dic = {}
for drugname in drugbank_dic['Generic_Name']:
    index = drugbank_dic['Generic_Name'].index(drugname)
    dbid = drugbank_dic['Primary_Accession_No'][index]
    drugname2dbid_dic[drugname] = dbid


# Extract the standard drug names (including synonyms) and cids:
with open(drugs_filename) as drugfile:
    name2cid_stdname_dic = {}
    lines = drugfile.readlines()
    for line in lines[1:len(lines)]:
        split_line = line.split('\n')[0].split('\t')
        std_name = split_line[0]
        synonyms = split_line[2].split('; ')
        pchem_cid = split_line[1]
        chembl_id = split_line[3]
        chebi_id = split_line[4]
        db_id = split_line[5]

        all_names = [std_name]+synonyms

        for name in all_names:
            name2cid_stdname_dic[name] = {}
            name2cid_stdname_dic[name]['cid'] = pchem_cid
            name2cid_stdname_dic[name]['name'] = std_name
            name2cid_stdname_dic[name]['chembl'] = chembl_id
            name2cid_stdname_dic[name]['chebi'] = chebi_id
            name2cid_stdname_dic[name]['db'] = db_id

# Open the file of interactions:
with open(infilename) as infile, open(outfilename, 'w') as outfile:
    writeHeader2File(outfile)
    lines = infile.readlines()
    for line in lines[1:len(lines)]:
        split_line = line.split('\n')[0].split('\t')
        for i in range(0, len(fields)):
            vars()[fields[i]] = split_line[i]

        # get all different identifiers, comma separated: pubchem_cid, chembl, chebi, drugbank
        pchem_cid = 'CID:'+ name2cid_stdname_dic[source]['cid']
        chembl_id = name2cid_stdname_dic[source]['chembl']
        chebi_id = name2cid_stdname_dic[source]['chebi']
        db_id = name2cid_stdname_dic[source]['db']
        if db_id == '': # for drugs retrieved from drugbank...
            try:
                db_id = drugname2dbid_dic[source]

            except:
                db_id = ''

        source_ids = ','.join(filter(None, [pchem_cid]+[chembl_id]+[chebi_id]+[db_id]))
        source_ids = re.sub(' ','',source_ids) # eliminate spaces if present
        compound_name = name2cid_stdname_dic[source]['name']
        if source != compound_name:
            print('\t' + source + ' to ' + compound_name)

        if re.match('"', compound_name):# eliminate quotes if present
            compound_name = compound_name[1:len(compound_name)-1]
        if re.match('"', source):
            source = source[1:len(source)-1]
        source_name = source #source name contains the original name from the article (in most cases = to compound name)
        writeLine2File(outfile)



