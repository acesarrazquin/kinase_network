#!/cm/shared/apps/python/2.7.6/bin/python

################################################################################################
#   Script to retrieve KSI from Newman et al. (PhosphoNetworks)
#   
#   Input: - comKSI.csv (kinase-substrate interactions)
#          - highResolutionNetwork (kinase-subs with phospho positions)
#
#   Output: phosphonetworks_ksi.txt -> has the KSI from comKSI with the phosphopositions in
#                the other file. Two more fields: "score" from comKSI, and "phospho_score" from the
#                other file.
#                
###############################################################################################
import argparse
import re
import mappings as map

# parse arguments
parser = argparse.ArgumentParser(description="Assembly of KSI from PhosphoNetworks (Newman et al.).")

parser.add_argument("phosphonetfile")
parser.add_argument("phosphonetfile2")
args = parser.parse_args()

# create dictionary of uniprot accession codes to gene symbols
uniprot2symbol, syn2uniprot = map.getUniprotMapDicts()

#create list of kinases uniprot accession codes
kinacc2name = map.getKinaseAcc2Symbol()

up_regex = "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"

PMIDs = '23549483'
dates = '2014'
sources = 'phosphonetworks'
outfields = ["source", "target", "source_name", "target_name", "PMIDs", "dates", "sources", "type", "phospho_positions"]
outfields.extend(["score", "phospho_score"])

notmapped = []
notkin = []
ksi_dic = {}
print("\nFormatting file %s...\n"%(args.phosphonetfile))
with open(args.phosphonetfile) as comksifile:
    lines = comksifile.readlines()
    for line in lines[1:len(lines)]:
        split_line = line.rstrip().split("\t")
        
        source_symb = split_line[0]
        target_symb = split_line[1]
        score = split_line[2]
        
        try:
            source =  syn2uniprot[source_symb]
        except:
            notmapped.append(source_symb)
            continue
        
        if not kinacc2name.has_key(source):
            notkin.append(source)
            continue

        source_name = uniprot2symbol[source]

        try:
            target = syn2uniprot[target_symb]
        except:
            notmapped.append(target_symb)
            continue

        target_name = uniprot2symbol[target]
        
        key = source + "||" + target
        if not ksi_dic.has_key(key):

            ksi_dic[key] = {}
        
            ksi_dic[key]["source"] = source
            ksi_dic[key]["target"] = target
            ksi_dic[key]["source_name"] = source_name
            ksi_dic[key]["target_name"] = target_name
            ksi_dic[key]["phospho_positions"] = []
            ksi_dic[key]["score"] = score
            ksi_dic[key]["PMIDs"] = PMIDs
            ksi_dic[key]["dates"] = dates
            ksi_dic[key]["sources"] = sources
            ksi_dic[key]["type"] = "KSI"
            ksi_dic[key]["phospho_score"] = []
        
        else:
            ksi_dic[key]["score"] = max(ksi_dic[key]["score"], score)

print("\nFormatting file %s...\n"%(args.phosphonetfile2))
with open(args.phosphonetfile2) as netfile:
    content = netfile.read().split(">")
    content = content[1:len(content)]

    for sub in content:
        lines = sub.rstrip().split("\n")

        target_symb = lines[0]
         
        try:
            target = syn2uniprot[target_symb]
        except:
            notmapped.append(target_symb)
            continue
        target_name = uniprot2symbol[target]

        for kin in lines[1:len(lines)]:
            split_kin = kin.split("\t") 

            source_symbol = split_kin[2]
            phospho_score = split_kin[3]
            phospho_positions = split_kin[1][1:len(split_kin[1])]

            # get uniprot accessions

            try:
                source = syn2uniprot[source_symbol]
            except:
                notmapped.append(source_symbol)
                continue
            
            if not kinacc2name.has_key(source):
                notkin.append(source)
                continue
            
            source_name = uniprot2symbol[source]
            
            key = source + "||" + target
            if not ksi_dic.has_key(key):
                print("NEw!!")

                ksi_dic[key] = {}
                
                ksi_dic[key]["source"] = source
                ksi_dic[key]["target"] = target
                ksi_dic[key]["source_name"] = source_name
                ksi_dic[key]["target_name"] = target_name
                ksi_dic[key]["phospho_positions"] = [phospho_positions]
                ksi_dic[key]["score"] = phospho_score
                ksi_dic[key]["PMIDs"] = PMIDs
                ksi_dic[key]["dates"] = dates
                ksi_dic[key]["sources"] = sources
                ksi_dic[key]["type"] = "KSI"
                ksi_dic[key]["phospho_score"] = phospho_score
            else:
                ksi_dic[key]["phospho_positions"].append(phospho_positions)
                ksi_dic[key]["phospho_score"].append(phospho_score)


notkin = list(set(notkin))
notmapped = list(set(notmapped))

with open("phosphonetworks_ksi.txt", "w") as outfile:
    
    for field in outfields:
        outfile.write(field + "\t")
    outfile.write("\n")
    
    for key in ksi_dic.keys():
        for field in outfields:
            value = ksi_dic[key][field]
            if type(value) == list:
                outfile.write(','.join(value) + "\t")
            else:
                outfile.write(value + "\t")
        outfile.write('\n')

print("\nSUMMARY:\n")
print("\tNot mapped by name to Uniprot: %d being %s\n" %(len(notmapped),','.join(notmapped)))
print("\tKinases not in list: %d being %s\n"%(len(notkin), ','.join(notkin)))
