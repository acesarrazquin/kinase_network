#!/cm/shared/apps/python/2.7.6/bin/python

###################################################################################################
# Script to retrieve ksi from PhosphoSitePlus.
#
# Input: -kinase_substrate_dataset, from Downloads
# Output: - file in database format with one more column "method" (invitro, invivo, both)
#
# Very easy, everything was using uniprot identifiers...
##################################################################################################

import argparse
import re
import mappings as map


# parse arguments
parser = argparse.ArgumentParser(description="Assembly of KSI from PhosphoSitePlus.")

parser.add_argument("phosphositefile")

args = parser.parse_args()

# create dictionary of uniprot accession codes to gene symbols
uniprot2symbol, syn2uniprot = map.getUniprotMapDicts()

#create list of kinases uniprot accession codes
kinacc2name = map.getKinaseAcc2Symbol()

up_regex = "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"

# read file
ksi_dic = {}
notmapped_up = []
notkin = []
with open(args.phosphositefile) as phosfile:
    lines = phosfile.readlines()
    for line in lines[1:len(lines)]:
        split_line = line.split("\t")
        # filter for human kinase and substrate
        if split_line[4] != "human" or split_line[10] != "human":
            continue

        # substrate
        subs_acc = split_line[7].split("-")[0]
        subs_name = split_line[5]
        subs_symb = split_line[8]
        subs_entrez = split_line[6]
        
        try:
            target_name = uniprot2symbol[subs_acc]
        except:
            try:
                subs_acc = syn2uniprot[subs_symb]
            except:
                notmapped_up.append(subs_acc)
                continue
            target_name = uniprot2symbol[subs_acc]
        target = subs_acc

        # kinase
        kin_acc = split_line[1].split("-")[0]
        kin_name = split_line[0]
        kin_symb = split_line[2]
        
        try:
            source_name = uniprot2symbol[kin_acc]
        except:
            try:
                kin_acc = syn2uniprot[kin_symb]
            except:
                notmapped_up.append(kin_acc)
                continue

            if not kinacc2name.has_key(kin_acc):
                notkin.append(kin_acc)
                continue
            
            source_name = uniprot2symbol[kin_acc]
        source = kin_acc
        
        # rest
        phospho_position = split_line[11][1:len(split_line[11])]
        
        if split_line[14] == 'X':
            if split_line[15] == 'X':
                method = "both"
            else:
                method = "in-vivo"
        elif split_line[15] == 'X':
            method = "in-vitro"
        else:
            method = ""

        # dictionary:
        key = source + "||" + target
        if not ksi_dic.has_key(key):
            ksi_dic[key] = {}
            
            ksi_dic[key]['source'] = source
            ksi_dic[key]['target'] = target
            ksi_dic[key]['source_name'] = source_name
            ksi_dic[key]['target_name'] = target_name
            ksi_dic[key]['PMIDs'] = "22135298"
            ksi_dic[key]['phospho_positions'] = [phospho_position]
            ksi_dic[key]['method'] = [method]
            ksi_dic[key]['dates'] = '2014'
            ksi_dic[key]['sources'] = "phosphositeplus"
            ksi_dic[key]['type'] = "KSI"
        else:
            ksi_dic[key]['phospho_positions'].append(phospho_position)
            ksi_dic[key]['method'].append(method)

notkin = list(set(notkin))
notmapped_up = list(set(notmapped_up))

# write file

outfields = ["source", "target", "source_name", "target_name", "PMIDs", "dates", "sources", "type", "phospho_positions"]
outfields.append("method")

print("\nWriting output file...\n")

with open("phosphositeplus_ksi.txt", "w") as outfile:
        for field in outfields:
            outfile.write(field + "\t")
        outfile.write("\n")

        for key in ksi_dic.keys():
            for field in outfields: 
                value = ksi_dic[key][field]
                if type(value) == list:
                    value = ",".join(set(value))
                outfile.write(value + "\t")
            outfile.write("\n")

print("\nSUMMARY:\n")
print("\tNot mapped to name by Uniprot: %d being %s\n" %(len(notmapped_up),','.join(notmapped_up)))
print("\tKinases not in list: %d being %s\n"%(len(notkin), ','.join(notkin)))
