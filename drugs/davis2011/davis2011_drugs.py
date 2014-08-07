#!/usr/bin/python

####################################################################################################
# File to extract the drug2protein interactions from Davis2011 paper.
#	- protein names/refseq ID are mapped to UniProtKB IDs, only for proteins in the initial
#	  kinase list
#	- drug names are mapped to PubChem Compound IDs using  Entrez utils from BioPython to query 
#	  PubChem with the alternative name, if exists, otherwise with the generic name
#     If several CIDs are retrieved, all of them will be considered (need manual curation)
#	- interactions are filtered to include only those with Kd <= 1000nM 
#
# Input: automatic (list of kinases, davis2011 compounds, davis2013 interactions and DrugBank files)
# 
# Output: file for each protein and drug interaction per line, tab separated  
#
####################################################################################################


def writeLine2File(outfile):
    outfile.write(source+'\t'+uniprot_id+'\t'+source_name+'\t'+prot_name+'\t'+
                  pmid+'\t'+dates+'\t'+type_+'\t'+pchem_cid+'\t'+target_form+'\t'+
                  strength+'\t'+strength_type+'\t'+strength_direction+'\t'+strong_target+'\t'+
                  technology+'\t'+cell_type+'\t'+sources+'\n')

def writeHeader2File(outfile):
    outfile.write('source'+'\t'+'target'+'\t'+'source_name'+'\t'+'target_name'+'\t'+'PMIDs'+'\t'+
                  'dates'+'\t'+'type'+'\t'+'source_id'+'\t'+'target_form'+'\t'+
                  'strength'+'\t'+'strength_type'+'\t'+'strength_direction'+'\t'+'strong_target'+'\t'+
                  'technology'+'\t'+'cell_type'+'\t'+'sources'+'\n')

#//////////////////////////////////// MAIN /////////////////////////////////////////////////////////
from drugs import *
import datetime

today = datetime.datetime.now().strftime('%Y%m%d')
outfilename = 'davis2011_drugs_formatted-'+today+'.txt'


compound_names_dic = {} # map compound name to alternative name
with open('davis2011_compounds.txt') as compounds_file:
    lines = compounds_file.readlines()
    for line in lines[1:len(lines)]:
        split_line = line.rstrip().split('\t')
        compound_name = split_line[0]
        alternative_name = split_line[1]

        if compound_name == 'INCB018424':
            compound_name = 'INCB18424' # this one is differently written in each file :S
        compound_names_dic[compound_name] = alternative_name

print('\nGetting UniProtID and DrugBank mappings...')
# mapping dictionaries
name2uniprot_dic = getName2UniprotDic()
uniprot2name_dic = getUniprot2NameDic()
ref2uniprot_dic = getRefseq2UniprotDic(name2uniprot_dic)
drugBank_dic = getDrugBankDic()

# constant field for output file:
strength_type = 'Kd[nM]'
strength_direction = 'low'
pmid = '22037378'
sources = pmid
dates = '2011'
type_ = 'DPI' # ?
technology = 'competition binding assay'
cell_type = ''

compound_targets_list = []
with open('davis2011_interactions.txt') as interactions_file, open(outfilename, 'w') as outfile:
    writeHeader2File(outfile)

    lines = interactions_file.readlines()

    header_fields = lines[0].rstrip().split('\t')
    compounds = header_fields[3:len(header_fields)]

    # map the compound names to PChemCID (this could be another identifier)
    print('\nMapping compound names to PubChem Compound identifiers...')
    name2pchemcid_dic = {}
    for compound_name in compounds: # Map using Entrez, if possible with the alternative name, otherwise the general name
        if compound_names_dic[compound_name]:
            alternative_name = compound_names_dic[compound_name]
            cids = getCIDFromName(alternative_name)
        else:
            cids = getCIDFromName(compound_name)
        pchem_cid = '; '.join(cids)
        name2pchemcid_dic[compound_name] = pchem_cid

    if pchem_cid == '':
        print('\t' + compound_name + ': no PubChem CID found!')

        # parse each protein and drug interactions for it in the file
    print('\nParsing drug-kinase interactions and writing output...')
    for line in lines[1:len(lines)]:
        split_line = line.rstrip().split('\t')

        ref_id = split_line[0]
        gene_name = split_line[1]
        kin_name = split_line[2] # this name contains the modifications in the protein (phosphorilation, mutation...)
        # map uniprot identifiers (based on the initial kinase list -> kinases not there are not mapped to Uniprot ID)
        try:
            uniprot_id = name2uniprot_dic[gene_name]
        except:
            try:
                uniprot_id = ref2uniprot_dic[ref_id]
                print('\t' + gene_name+': UniProt_ID recovered by RefSeq ID' )
            except:
                if re.match('[PQ]', ref_id):
                    uniprot_id = ref_id.split('.')[0]
                else:
                    uniprot_id = ''
                    print('\t' + gene_name + ' (' + ref_id + ') : not mapped to UniProtID')
                continue
        if uniprot_id != '':
            prot_name = uniprot2name_dic[uniprot_id]
        else:
            prot_name = gene_name

        if re.search('\(', kin_name):
            target_form = re.sub('\)', '', kin_name.split('(')[1])
        elif re.search('-nonph', kin_name) or \
                re.search('-ph', kin_name) or \
                re.search('-cyc', kin_name):
            target_form = kin_name.split('-')[1]
        else:
            target_form = ''

        if re.search('"', target_form):
            target_form = target_form[0:len(target_form)-1]


        kd_values = split_line[3:len(split_line)]
        for i in range(0, len(kd_values)):
            kd_val = kd_values[i]
            if kd_val and float(kd_val) <= 1000: # filter out interactions with kd > 1000nM
                compound_name = compounds[i]
                try:
                    pchem_cid = name2pchemcid_dic[compound_name]
                except:
                    pchem_cid = ''

                if float(kd_val) <= 100:
                    strong_target = 'yes'
                else:
                    strong_target = 'no'

                strength = kd_val

                if compound_names_dic[compound_name]:
                    source = compound_names_dic[compound_name]
                else:
                    source = compound_name
                source_name = compound_name
                writeLine2File(outfile)
