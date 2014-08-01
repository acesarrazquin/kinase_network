#!/cm/shared/apps/python/2.7.6/bin/python

#########################################################################################################################3
# Script to extract PPI from Reactome file "Human protein-protein interaction pairs in tab-delimited format"
#
# Input: - reactome file
#
# Output: -reactome_ppi.txt, in database format with some added columns (reaction_type, reaction_id) 
#
# Possibitlity to filter for interaction type:
#       - association
#       - physical association
########################################################################################################################
import argparse
import re
from xml.etree import ElementTree as et
import mappings as map

# parse arguments
parser = argparse.ArgumentParser(description="Extraction of KSI from Reactome interaction file, in mitab format.")

parser.add_argument("reactomefile")
parser.add_argument("-t", "--int_type", help="interaction type", choices=["association", "physical association", "any"], default="any")

args = parser.parse_args()

# create list of kinases uniprot accession codes
uniprot2symbol, syn2uniprot = map.getUniprotMapDicts()

# read reactome interactions file
up_regex = "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
ppi_dic = {}
print("\nReading Reactome Interactions file...\n")
with open(args.reactomefile) as rfile:
    lines = rfile.readlines()
    for line in lines[1:len(lines)]:
        
        split_line = line.rstrip().split("\t")

        source = split_line[0].split(":")[1].split("-")[0]
        target = split_line[1].split(":")[1].split("-")[0]
        
        # only uniprot proteins
        if not re.match(up_regex, source) or not re.match(up_regex, target):
            continue

        # apply interaction type filter
        int_type = split_line[11].split("(")[1]
        int_type = int_type[0:len(int_type)-1]

        if args.int_type != "any" and args.int_type != int_type:
            continue

        # no self-interactions # for the moment I include everything
        #if source==target:
        #    continue

        complex_ids = split_line[17].split("|")
        #complex_ids = []
        #for complex in complex_ids:
        #    complex_id = complex.split(":")[1]
        #    complex_ids.append(complex_id)

        pmid = split_line[8].split(":")[1] # there is only the general for reactome...

        key = source + "||" + target
        if key not in ppi_dic:
            ppi_dic[key] = {}
            
            ppi_dic[key]["source"] = source
            ppi_dic[key]["target"] = target
            ppi_dic[key]["interaction_type"] = [int_type]
            ppi_dic[key]["complex_id"] = complex_ids
            ppi_dic[key]["pmid"] = [pmid]
        else:
            ppi_dic[key]["interaction_type"].extend([int_type])
            ppi_dic[key]["complex_id"].extend(complex_ids)
            ppi_dic[key]["pmid"].extend([pmid])

outfilename = "reactome_mitab_ppi_" + args.int_type + ".txt"
outfilename = re.sub(" ", "", outfilename)
fields = ["source", "target", "source_name", "target_name", "PMIDs", "dates", "sources", "type", 
          "source_is_bait", "target_is_bait"]
fields.extend(["complex_id", "interaction_type"]) # in case i want to add them
dates = '2014'
sources = 'reactome'
type = 'PPI'
source_is_bait = "no"
target_is_bait = "no"
nomap = []
print("\nWriting output file...\n")
with open(outfilename, 'w') as outfile:
    
    for field in fields:
        outfile.write(field + "\t")
    outfile.write("\n")

    for key in ppi_dic.keys():
        source = ppi_dic[key]["source"]
        target = ppi_dic[key]["target"]

        try:
            source_name = uniprot2symbol[source]
        except:
            nomap.append(source)
            continue

        try:
            target_name = uniprot2symbol[target]
        except:
            nomap.append(target)
            continue

        PMIDs = ','.join(set(ppi_dic[key]["pmid"]))

        interaction_type = ','.join(set(ppi_dic[key]["interaction_type"]))
        complex_id = ','.join(set(ppi_dic[key]["complex_id"]))

        for field in fields:
            outfile.write(vars()[field] + "\t")
        outfile.write("\n")

print("\nProteins not mapped to Uniprot: %d being %s"%(len(set(nomap)), ','.join(set(nomap))))
