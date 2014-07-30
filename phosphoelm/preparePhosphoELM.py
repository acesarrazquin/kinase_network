#!/cm/shared/apps/python/2.7.6/bin/python

#########################################################################################################
#   Script to retrieve KSI from PhosphoELM dataset.   
#   
#   Input: - dataset downloaded from website (I first extract only human (grep "Homo sapiens") from the vertebrates dataset
#          - kinase names files, also from website phophoelm.../kinases.html, to map kinase names to uniprot
#   Output: - file in database format with an extra field "experiment" (LTP/HTP, low or high throughput)
#   
#   Notes: - needed to
#             -convert ensp to uniprot (because some substrates are written as ensp)
#             -convert uniprot to symbol
#             -secondary accession to primary accession? todo
#          - for many substrates there is no kinase information
#       
#########################################################################################################
import argparse
import re
import mappings as map

# parse arguments
parser = argparse.ArgumentParser(description="Assembly of KSI from PhosphoELM.")

parser.add_argument("phosphoelmfile")
parser.add_argument("kinasenamefile")

args = parser.parse_args()


# create dictionary of uniprot accession codes to gene symbols
uniprot2symbol, syn2uniprot = map.getUniprotMapDicts()

# create a dictionary of ensp identifiers to uniprot accession codes
ensp2uniprot = map.getEnsp2Uniprot()

#create list of kinases uniprot accession codes
kinacc2name = map.getKinaseAcc2Symbol()

up_regex = "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"

# read kinase names file:
print("\nReading kinase names file...\n")
phosphoelmkin2up = {}
notfound_kin = []
with open(args.kinasenamefile) as kinfile:
    lines = kinfile.readlines()
    for line in lines[1:len(lines)]:
        split_line = line.split("\t")

        kinase_up = split_line[2].split("-")[0]
        kinase_name = split_line[0]

        if re.search('group', kinase_name):
            continue
        elif kinase_up == "":
            # try to match by name
            try:
                kinase_up = syn2uniprot[kinase_name]
            except:
                notfound_kin.append(kinase_name)
        else:
            phosphoelmkin2up[kinase_name] = kinase_up
notfound_kin = list(set(notfound_kin))

# read phosphoELM file:
print("\nReading PhosphoELM file...\n")

ksi_dic = {}

multimapped_ensp = []
notmapped_ensp = []
notmapped_up = []
nolist_kin = []
with open(args.phosphoelmfile) as elmfile:
    lines = elmfile.readlines()
    for line in lines[0:len(lines)]:
        split_line = line.split("\t")
        
        if split_line[7] != "Homo sapiens":
            continue

        # substrate
        substrate_acc = split_line[0].split("-")[0] # can be uniprot or ensp
        
        if re.search("ENSP0", substrate_acc):
            try:
                substrate_up = ensp2uniprot[substrate_acc]
            except:
                notmapped_ensp.append(substrate_acc)
                continue
            if len(substrate_up) > 1:
                multimapped_ensp.append(substrate_acc)
                continue
            else:
                substrate_acc = substrate_up[0]
        elif not re.match(up_regex, substrate_acc):
            continue

        try:
            target_name = uniprot2symbol[substrate_acc]
        except: 
            #print("UniProt not in list: %s"%(substrate_acc)) # not mapped
            notmapped_up.append(substrate_acc)
            continue
        target = substrate_acc
        
        # kinase
        kinase_name = split_line[5]
        
        try:
            source = phosphoelmkin2up[kinase_name]
        except:
            continue

        try:
            source_name = kinacc2name[source]
        except:
            nolist_kin.append(source)
            continue
        
        # rest of fields
        phospho_position = split_line[2] # only one in each line
        pmid = split_line[4] # only one in each line, can be N.N.
        if pmid == "N.N.":
            pmid == ""
        experiment = split_line[6] # only one
        dates = split_line[8][0:4]

        # dictionary:
        key = source + "||" + target
        if not ksi_dic.has_key(key):
            ksi_dic[key] = {}

            ksi_dic[key]['source'] = source
            ksi_dic[key]['target'] = target
            ksi_dic[key]['source_name'] = source_name
            ksi_dic[key]['target_name'] = target_name
            
            ksi_dic[key]['PMIDs'] = [pmid]
            ksi_dic[key]['phospho_positions'] = [phospho_position]
            ksi_dic[key]['experiment'] = [experiment]
            ksi_dic[key]['dates'] = [dates]
            ksi_dic[key]['sources'] = "phospho.elm"
            ksi_dic[key]['type'] = "KSI"
        else:
            ksi_dic[key]['PMIDs'].append(pmid)
            ksi_dic[key]['phospho_positions'].append(phospho_position)
            ksi_dic[key]['experiment'].append(experiment)
            ksi_dic[key]['dates'].append(dates)

multimapped_ensp = list(set(multimapped_ensp))
notmapped_ensp = list(set(notmapped_ensp))
notmapped_up = list(set(notmapped_up))
nolist_kin = list(set(nolist_kin))

# write file
print("\nWriting output file...\n")

outfields = ["source", "target", "source_name", "target_name", "PMIDs", "dates", "sources", "type", "phospho_positions"]
outfields.append("experiment")

with open("phosphoelm_ksi.txt", 'w') as outfile:
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
print("\tNot mapped from ENSP to UniProt: %d being %s\n" %(len(notmapped_ensp), ','.join(notmapped_ensp)))
print("\tMultimapped from ENSP to UniProt: %d being %s\n"%(len(multimapped_ensp), ','.join(multimapped_ensp)))
print("\tNot mapped to name by Uniprot: %d being %s\n" %(len(notmapped_up),','.join(notmapped_up)))
print("\tKinases not in list: %d being %s\n"%(len(nolist_kin), ','.join(nolist_kin)))
