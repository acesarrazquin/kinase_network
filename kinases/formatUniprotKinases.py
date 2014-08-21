#!/cm/shared/apps/python/2.7.6/bin/python

###########################################################################################################################
# To format the file with kinases retrieved from uniprot: uniprot.org/docs/pkinfam.txt
#
#
#

import re
import argparse
import mappings as map
import datetime
import copy

# parse arguments
parser = argparse.ArgumentParser(description="Formatting uniprot human and mouse kinases file.")

parser.add_argument("kinasesfile")
parser.add_argument("-o", "--outputfile", help="name of output file")

args = parser.parse_args()

# get dictionaries
uniprot2symbol, syn2uniprot, uniprot2syn = map.getUniprotMapDicts(synonyms=True)
englishwords = map.getEnglishWords()

# read file
kinases_types = {} #dictionary for types and descriptions
kinases_list = []

atypical = 0
#type_pattern = '\*'
group_pattern = '===='
protein_pattern = '_HUMAN'

with open(args.kinasesfile) as kinfile:
    all_lines = kinfile.read()
    split_lines = all_lines.split('\r')

    for i in range(0, len(split_lines)):
        line = split_lines[i]

        if re.search(group_pattern, line) and re.search(group_pattern, split_lines[i+2]):
            kinase_type = re.split('\\n| protein| kinase| family| Ser/Thr', split_lines[i+1])[1]
            print(kinase_type)

            for j in range(i+3, len(split_lines)):
                line2 = split_lines[j]

                if re.search(protein_pattern, line2):
                    split_kinase = line2.split()
                    kinase_name = split_kinase[0]
                    kinase_uniprotid = split_kinase[2][1:len(split_kinase[2])] # care with this, in new file its different

                    synonyms = uniprot2syn[kinase_uniprotid]

                    # apply words filter:
                    filt_synonyms = copy.copy(synonyms)

                    for synonym in synonyms:
                        synonym_ = synonym.upper()

                        if synonym_ in englishwords:
                            # print(kinase_name + ' : ' + synonym + ' removed -- English word')
                            filt_synonyms.remove(synonym)

                    # add the initial name if not present
                    if kinase_name not in synonyms:
                        synonyms.insert(0, kinase_name)

                    synonyms = ', '.join(synonyms)
                    filtsynonyms = ', '.join(filt_synonyms)

                    # append to list
                    kinase = [kinase_uniprotid, kinase_name, synonyms, filtsynonyms, kinase_type]
                    kinases_list.append(kinase)

                if re.search(group_pattern, line2):
                    break


today = datetime.date.today().strftime("%Y%m%d")
if args.outputfile:
    outputfile = args.outputfile
else:
    outputfile = args.kinasesfile.split('.txt')[1]+'_formated'+today+'.txt'

with open(outputfile, 'w') as outfile:
    outfile.write('uniprot_acc\tname\tsynonyms\tfilt_synonyms\ttype\n')
    for line in kinases_list:
        outfile.write(line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+line[3]+'\t'+line[4]+'\n')

