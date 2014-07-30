#!/cm/shared/apps/python/2.7.6/bin/python

#########################################################################################################################3
# Script to extract PPI from Reactome file "Human protein-protein interaction pairs in tab-delimited format"
#
# Input: - reactome file
#        - genenamefile, to convert accession numbers to gene symbols
#
# Output: -reactome_ppi.txt, in database format with some added columns (reaction_type, reaction_id) 
#
# Possibitlity to filter for reaction type...Now really to the interaction types:
#       - direct_complex - interactors are directly in the same complex, i.e. 
#                           w/o further nested complexes. From the example interactions A <->B and 
#                           C<-D> are of this type.
#       - indirect_complex - interactors in different subcomplexes of a complex.
#                           From the example above interactions A<->C, A<->D, B<->C and B<->D are 
#                           of this type.
#       - reaction - interactors participate in the same reaction. Only those 
#                   reactions are reported for which the intreactors are not complexed (with 
#                   the exception being association dissociation reactions which are reported).
#       - neighbouring_reaction - interactor participate in 2 "consecutive"
#                               reactions, i.e. when one reaction produces a PhysicalEntity which is either an 
#                               input or a catalyst for another reaction. However, to avoid the "trivial" (i.e.
#                               over ATP, ADP etc) interactions, the computation is done using the
#                               'precedingEvent' attribute used to order reactions in a pathway.
########################################################################################################################
import argparse
import re
from xml.etree import ElementTree as et
import mappings as map

# parse arguments
parser = argparse.ArgumentParser(description="Extraction of KSI from Reactome SBML files.")

parser.add_argument("reactomefile")
#parser.add_argument("-o", "--outputfile", help="name of output file")
#parser.add_argument("-g", "--genenamefile", default="inputfiles/sp2name_human-20140703.tab", dest="gnfile",
 #                   help="table of uniprot accesion codes, gene symbols, protein names and descriptions.")
parser.add_argument("-i", "--interactiontype", choices=["direct", "complexes", "reaction", "all"], default="all",
                    help="type of interactions to consider (each choices includes previous)", dest="int_type")
args = parser.parse_args()

# create list of kinases uniprot accession codes
uniprot2symbol, syn2uniprot = map.getUniprotMapDicts()

# read reactome interactions file
up_regex = "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
ppi_dic = {}
print("\nReading Reactome Interactions file...\n")
with open(args.reactomefile) as rfile:
    for line in rfile:
        if line[0] == "#":
            continue

        
        split_line = line.rstrip().split("\t")
        source = split_line[0].split(":")[1].split("-")[0]
        target = split_line[3].split(":")[1].split("-")[0]
        
        # only uniprot proteins
        if not re.match(up_regex, source) or not re.match(up_regex, target):
            continue
        
        # no self-interactions
        if source==target:
            continue

        react_type = split_line[6]
        react_id = split_line[7]

        # filter interactions
        if args.int_type == "direct":
            if react_type != "direct_complex":
                continue
        if args.int_type == "complex":
            if react_type != "direct_complex" and react_type != "indirect_complex":
                continue
        if args.int_type == "reaction":
            if react_type == "neighbouring_reaction":
                continue

        try:
            pmids = split_line[8].split(",")
        except:
            pmids = []

        key = source + "||" + target
        if key not in ppi_dic:
            ppi_dic[key] = {}
            
            ppi_dic[key]["source"] = source
            ppi_dic[key]["target"] = target
            ppi_dic[key]["reaction_type"] = [react_type]
            ppi_dic[key]["reaction_id"] = [react_id]
            ppi_dic[key]["pmids"] = pmids
        else:
            ppi_dic[key]["reaction_type"].extend([react_type])
            ppi_dic[key]["reaction_id"].extend([react_id])
            ppi_dic[key]["pmids"].extend(pmids)

outfilename = "reactome_ppi.txt"
fields = ["source", "target", "source_name", "target_name", "PMIDs", "dates", "sources", "type", 
          "source_is_bait", "target_is_bait"]
fields.extend(["reaction_id", "reaction_type"]) # in case i want to add them
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

        PMIDs = ','.join(set(ppi_dic[key]["pmids"]))

        reaction_type = ','.join(set(ppi_dic[key]["reaction_type"]))
        reaction_id = ','.join(set(ppi_dic[key]["reaction_id"]))

        for field in fields:
            outfile.write(vars()[field] + "\t")
        outfile.write("\n")

print("\nProteins not mapped to Uniprot: %d being %s"%(len(set(nomap)), ','.join(set(nomap))))
