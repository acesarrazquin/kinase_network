#!/cm/shared/apps/python/2.7.6/bin/python

##############################################################################################################
# prepareNetworKIN.py changes the predictions made by NetworKIN v3 to the common database format.
#
# Necessary to map ENSP identifiers to UniProt and GeneNames
##############################################################################################################

import argparse
import re
import mappings as map

# parse arguments
parser = argparse.ArgumentParser(description="Assembly of KSI from NetworKIN predictions.")

parser.add_argument("networkinfile")
parser.add_argument("-o", "--outputfile", help="name of output file")

args = parser.parse_args()


# create dictionary of uniprot accession codes to gene symbols
uniprot2symbol, syn2uniprot = map.getUniprotMapDicts()

# create a dictionary of ensp identifiers to uniprot accession codes
ensp2uniprot = map.getEnsp2Uniprot()

# create list of kinases uniprot accession codes
kinacc2name = map.getKinaseAcc2Symbol()

# format networkin data
if args.outputfile:
    outputfile = args.outputfile
else:
    outputfile = args.networkinfile.split(".")[0]+"_formated.txt"

outfields = ["source", "target", "source_name", "target_name", "PMIDs", "dates", "sources", "type", "phospho_positions"]
outfields.extend(["networkin_score", "netphorest_score", "string_score"])
nomap_uniprot = []
kinnot = []
nomap_ensp = []
with open(args.networkinfile) as netfile, open(outputfile, "w") as outfile:
    print("\nFormatting NetworKIN data...\n")

    # write headers
    for field in outfields:
        outfile.write(field + "\t")
    outfile.write("\n")

    lines = netfile.readlines()
    counter = 0
    for line in lines:
        counter += 1
        if counter == 1:
            fields = re.sub("#", "", line).rstrip().split("\t")
            #outfile.write("\t".join(outfields)+"\n") # could use this strategy for all cases

        else:
            split_line = line.rstrip().split("\t")
            
            for i in range(0, len(fields)):
                vars()[fields[i]] = split_line[i]
            
            # get sources and targets (a bit messy in networKIN...)
            target_symbol = substrate_name
            target_ensp = string_identifier
            target_ensp_alt = re.sub("\)", "", substrate.split(" (")[1]) 
            
            source_symbol = id
            source_ensp = string_path.split(",")[0]
            
            # map ensp to uniprot and from there to gene symbols
            try:
                target = ensp2uniprot[target_ensp]
            except:
                try:
                    target = ensp2uniprot[target_ensp_alt] #try with the other field
                except:
                    try:
                        target = [syn2uniprot[target_symbol.upper()]] # try by name
                    except:
                        nomap_ensp.append(target_ensp)
                        continue    
            
            target_filt = [] # targets only if uniprot in list
            target_names = []
            for tar in target: # now it's a list
                try:
                    target_name = uniprot2symbol[tar]
                    target_names.append(target_name)
                    target_filt.append(tar)
                except:
                    nomap_uniprot.append(tar)
                    continue
            
            try:
                source = ensp2uniprot[source_ensp]
            except:
                try:
                    source = [syn2uniprot[source_symbol.upper()]]
                except:
                    nomap_ensp.append(source_ensp)
                    continue
           
            source_filt = []
            source_names = []
            for src in source:
                try:
                    source_name = uniprot2symbol[src]
                    source_names.append(source_name)
                    source_filt.append(src)
                except:
                    nomap_uniprot.append(src)
                    continue

                try:
                    kinacc2name[src] # to check is kinase
                except:
                    kinnot.append(src)
                    continue

                source_filt.append(src)
                source_names.append(source_name)
                    #try:
                        #if source_symbol=="PKCalpha":
                        #    source_symbol = "PRKCA"
                        #elif source_symbol=="CK1alpha":
                        #    source_symbol = "CSNK1A1"
                        #elif source_symbol=="IKKbeta":   
                        #    source_symbol = "IKBKB"                                                                                                                                             
                    #    source = syn2uniprot[source_symbol.upper()]
                    #    source_name = kinacc2name[src]
                    #    source_names.append(source_name)
                    #except:
                    
                   #     kinnot.append(source_symbol)
                    #    continue

            for i in range(0, len(source_filt)):
                for j in range(0, len(target_filt)):

                    source = source_filt[i]
                    source_name = source_names[i]
                    target = target_filt[j]
                    target_name = target_names[j]
                
                    phospho_positions = position

                    # other fields:
                    PMIDs = "24874572"
                    sources = "networkin"
                    dates = "2014"
                    type = "KSI"

                    for field in outfields:
                        outfile.write(vars()[field] + "\t")
                    outfile.write("\n")

print("\nSUMMARY:\n")
print("\tNot mapped from ENSP to UniProt: %d being %s" %(len(set(nomap_ensp)), ','.join(set(nomap_ensp))))
print("\tNot mapped to name by Uniprot: %d being %s" %(len(set(nomap_uniprot)),','.join(set(nomap_uniprot))))
print("\tKinases not in list: %d being %s"%(len(set(kinnot)), ','.join(set(kinnot))))

