#!/cm/shared/apps/python/2.7.6/bin/python
############################################################################################################################
# prepareInteractDB.py parses a file extracted from our database of protein-protein interactions (by Peter), which already
# applies some filters (e.g. no biogrid, only physical interactions, etc...), and generates an output file in an adequate format
# to generate the Human Protein Kinase Network.
#
# To extract PPI or KSI ("phosphorylation reaction"), with possibility to include CeMM confidential PPIs as well as from
# Couzens et al. paper (though they are already included in the latest IntAct version)
#
#
# Input: interactionfile containing fields [p1_accession, p2_accession, p1_database, p2_database, p1_name, p1_taxid, p2_name, p2_taxid,
#                                           p1_molecule_type, p2_molecule_type, experimentalrole, i_accs, sources, detectionmethods, 
#                                           interactiontypes, publications, nexp]   
#
# Output: file with fields ["source", "target", "source_name", "target_name", "PMIDs", "dates", "sources", "type"]
#   if PPI also ["source_is_bait", "target_is_bait"]
#   if KSI also ["phospho_positions"] (although this field will be empty for this dataset)
#
###########################################################################################################################
import argparse
import re
import mappings as map

# parse arguments
parser = argparse.ArgumentParser(description="Assembly of interactions from InteractionDB files.")

parser.add_argument("interactionfile")
parser.add_argument("-o", "--outputfile", help="name of output file")
parser.add_argument("-t", "--type", choices=["PPI", "KSI"], default="PPI", help="type of interaction data")
parser.add_argument("-x", "--taxon", choices=["human", "mouse"], default="human")
parser.add_argument("-c", "--withcemm", action="store_true", help="include CeMM confidential interactions")
parser.add_argument("-z", "--couzens", action="store_true", help="include interactions from 'Couzens et al.'")

args = parser.parse_args()

# change taxon name to ID
if args.taxon == "human":
    taxon_id = 9606
elif args.taxon == "mouse":
    taxon_id = 10090

# create dictionary of uniprot accession codes to gene symbols
uniprot2symbol, syn2uniprot = map.getUniprotMapDicts()


# create list of kinases uniprot accession codes
kinacc2name = map.getKinaseAcc2Symbol()

# create dictionary of interactions
notfound = []
interact = {}
with open(args.interactionfile) as intfile:
    print("\nRetrieving interactions from file...\n")

    counter = 0
    for line in intfile:
        counter += 1

        split_line = line.rstrip().split("\t")

        if counter == 1: 
            fields = split_line

        else:

            # GET VARIABLES and APPLY FILTERS
            
            ## get variables and values
            for i in range(0, len(fields)): 
                vars()[fields[i]] = split_line[i]
            
            interactiontypes_list = interactiontypes.split(",")
            sources_list = sources.split(",")
            publications_list = publications.split("@@")
            experimentalrole_list = experimentalrole.split(",")
            
            # KSI vs PPI (if interactiontype contains 'phosphorylation reaction', consider only as KSI)
            if "phosphorylation reaction" in interactiontypes_list:
                if args.type == "PPI":
                    continue
            else:
                if args.type == "KSI":
                    continue
            
            # UniProt pairs only
            if p1_database != "UniProt" or p2_database != "UniProt":
                continue


            # 1. do not distinguish between isoforms: "-[0-9]" at the end of accesssion code
            # 2. remove naming errors detected in some cases...
            isoform_pattern = "-[1-9]$"
            wrong_pattern1 = "PRO_"
            wrong_pattern2 = "\|$"
            
            p_from = re.sub(isoform_pattern,"", p1_accession)
            p_from = re.sub(wrong_pattern1,"", p_from)
            p_from = re.sub(wrong_pattern2,"", p_from)

            p_to = re.sub(isoform_pattern, "", p2_accession)
            p_to = re.sub(wrong_pattern1,"", p_to)
            p_to = re.sub(wrong_pattern2,"", p_to)


            # Maps Uniprot accessions to gene names -> if not, then continue (not in my Swissprot list)
            try:
                p_from_name = uniprot2symbol[p_from]
            except:
                notfound.append(p_from)
                continue
                #print("protein %s not found in dictionary"%(p_from))

            try:
                p_to_name = uniprot2symbol[p_to]
            except:
                notfound.append(p_to)
                continue
                #print("protein %s not found in dictionary"%(p_to))

            # No self interactions
            if p_from == p_to:
                continue

           # Check for unambiguous KSI.
               # NOTE: the "phosphorylation reaction" annotation comes from IntAct and it's ambiguous in the directionality.
            kinases = kinacc2name.keys()
            if args.type == "KSI":
                if (p_from not in kinases) and (p_to not in kinases): # inconsistency
                    print("\tNo kinase in KSI: %s %s"%(p_from, p_to))
                    continue
                if (p_from in kinases) and (p_to in kinases): # not possible to know which of the two is a substrate
                    print("\tBoth are kinases, ambiguous kinase-substrate interaction!: %s %s"%(p_from, p_to))
                    continue
                if (p_to in kinases): # if only one is a kinase it cannot be the substrate! reverse order:
                    print("\tWrong order in KSI: %s %s"%(p_from, p_to))
                    tmp = p_from
                    p_from = p_to
                    p_to = tmp
            
           # apply CeMM confidential interactions filter:
            if not args.withcemm:
                if "cemm_confidential" in sources_list:
                    if len(sources_list) == 1:
                        continue
                    else:
                        while 1:
                            try:
                                ind = sources_list.index("cemm_confidential")
                            except:
                                break
                            del sources_list[ind]
                            del publications_list[ind]

           # Get PMIDs and dates (years)
            pmids = []
            dates = [] # only the year
            for publication in publications_list:
                publication_split = publication.split("~")
                pmid = publication_split[0]
                pmid = re.sub("^MINT-", "", pmid)
                ## what does this do here? $elem =~ /^([0-9]+)~[A-Z]/;
                if pmid == "N/A" or pmid == "":
                    pmid = "NA"
                elif re.search("unassigned", pmid):
                    pmid = "NA"
                pmids.append(pmid)
                
                year = publication_split[2]
                if year == "N/A":
                    year = "NA"
                elif year != "in preparation":
                    year = year[0:4]
                dates.append(year)

            # CREATE DICTIONARY 
            
            # the dictionary key has both proteins sorted by name 
            sorted_p = sorted([p_from, p_to])
            key = sorted_p[0] + "||" + sorted_p[1]

            if key not in interact.keys(): # create key and fields
                interact[key] = {}

                interact[key]["source"] = p_from
                interact[key]["target"] = p_to
                interact[key]["source_name"] = p_from_name
                interact[key]["target_name"] = p_to_name
                interact[key]["source_is_bait"] = "no"
                interact[key]["target_is_bait"] = "no"
                interact[key]["sources"] = sources_list
                interact[key]["PMIDs"] = pmids
                interact[key]["dates"] = dates
                interact[key]["type"] = args.type
                interact[key]["phospho_positions"] = ""

            else: # append data of sources, pmids, and dates
                for i in range(0,len(sources_list)):
                    if (pmids[i] not in interact[key]["PMIDs"]) or (sources_list[i] not in interact[key]["sources"]):
                        interact[key]["sources"].append(sources_list[i])
                        interact[key]["PMIDs"].append(pmids[i])
                        interact[key]["dates"].append(dates[i])
                        #print("Data appended for key: %s"%(key))

            # Assigns roles checking possible reversed source/target in previous lines
            for role in experimentalrole_list:
                if role == "bait-prey":
                    if p_from == interact[key]["source"]:
                        interact[key]["source_is_bait"] = "yes"
                    else:
                        interact[key]["target_is_bait"] = "yes"
                elif role == "prey-bait":
                    if p_from == interact[key]["source"]:
                        interact[key]["target_is_bait"] = "yes"
                    else:
                        interact[key]["source_is_bait"] = "yes"


            
# data from Couzens et al. "Protein Interaction Network of the Mammalian Hippo Pathway Reveals Mechanisms of Kinase-Phosphatase Interactions"
# IT SHOULD BE INCLUDED IN NEW UPDATE OF THE PPI DATABASE!!
if args.couzens:
    with open("couzens-et-al.txt") as couzensfile:
        print("\nRetrieving interactions from Couzens et al...\n")
        
        counter=0
        for line in couzensfile:
            counter +=1 
            if counter == 1:
                continue
            
            split_line = line.split("\t")
            
            p_from = split_line[0].split(":")[1]
            p_to = split_line[1].split(":")[1]

            # Maps Uniprot accessions to gene names
            try:
                p_from_name = uniprot2symbol[p_from]
            except:
                notfound.append(p_from)
                continue

            try:
                p_to_name = uniprot2symbol[p_to]
            except:
                notfound.append(p_to)
                continue

            # no self interaction
            if p_from == p_to:
                continue
            
            # append to dictionary if new
            sorted_p = sorted([p_from, p_to])
            key = sorted_p[0] + "||" + sorted_p[1]

            if key not in interact.keys(): # create key and fields
                interact[key] = {}

                interact[key]["source"] = p_from
                interact[key]["target"] = p_to
                interact[key]["source_name"] = p_from_name
                interact[key]["target_name"] = p_to_name
                interact[key]["source_is_bait"] = "no"
                interact[key]["target_is_bait"] = "no"
                interact[key]["sources"] = "intact"
                interact[key]["PMIDs"] = "24255178"
                interact[key]["dates"] = "2013"
                interact[key]["type"] = args.type
            else:
                if ("intact" not in interact[key]["sources"]) or ("24255178" not in interact[key]["PMIDs"]):
                    interact[key]["sources"].append("intact")
                    interact[key]["PMIDs"].append("24255178")
                    interact[key]["dates"].append("2013")

            if split_line[18] == 'psi-mi:"MI:0496"(bait)':
                interact[key]["source_is_bait"] = 'yes'
                      
            if split_line[19] == 'psi-mi:"MI:0496"(bait)':
                interact[key]["target_is_bait"] = 'yes'



# output file:
if args.outputfile:
    outfilename = args.outputfile
else:
    outfilename = args.interactionfile.split(".")[0]
    outfilename += "_" + args.taxon + "_" + args.type
    if args.withcemm:
        outfilename += "withcemm"
    if args.couzens:
        outfilename += "couzens"
    outfilename += ".txt"

print("\nWriting output file %s...\n"%(outfilename))
with open(outfilename, "w") as outfile:
    
    fields = ["source", "target", "source_name", "target_name", "PMIDs", "dates", "sources", "type"]
    if args.type == "PPI":
        fields.extend(["source_is_bait", "target_is_bait"])
    elif args.type == "KSI":
        fields.append("phospho_positions")


    # write header
    for field in fields:
        if fields.index(field) == len(fields)-1:
            outfile.write(field+"\n")
        else:
            outfile.write(field+"\t")

    # write contents
    for key in interact.keys():
        for field in fields:
            towrite = interact[key][field]
            if fields.index(field) == len(fields)-1:
                if type(towrite) == list:
                    towrite = set(towrite)
                    outfile.write(",".join(towrite) + "\n")
                else:
                    outfile.write(towrite + "\n")
            else:
                if type(towrite) == list:
                    towrite = set(towrite)
                    outfile.write(",".join(towrite) + "\t")
                else:
                    outfile.write(towrite + "\t")

notfound = set(notfound)
print("\nNOTE: %d protein accession numbers not mapped to gene name: %s\n"%(len(notfound), ','.join(notfound)))











