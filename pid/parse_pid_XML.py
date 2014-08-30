#!/cm/shared/apps/python/2.7.6/bin/python
###################################################################################################################################
# Script to retrieve PPI and KSI from PID XML file "NCI-Nature_Curated.xml" and similar (from BioCarta, from Reactome)
# 
# Uses as input: - DIP XML file
#
# Gets first proteins, then families, then complexes (in that order, to follow the structure of the file)
# PPI are considered between all elements in the complex.
# 
# Checks if it's phosphorylation by comparing the modifications in inputs and outputs of "modification" interactions.
# Phosphorylation positions are included.
# Filters kinases so that we only include the ones in the kinase list. If several kinases are involved in the reaction, each
# KSI (one for each kinase) is considered, but are marked as "expanded" (possible to filter).
#
# Evidence codes also included: (FILTER FOR EVIDENCE???)
#   - IAE (Inferred from Array Experiments), IC (Inferred by Curator), IDA (Direct Assay), IFC (Funcional Complementation),
#     IGI (Genetic Interaction), IMP (Mutant Phenotype), IOS (Other Species), IPI (Physical Interaction), RCA (Reviewed COmputational Analysis),
#     RGE (Reporter Gene Expression), TAS (Traceable Author Statement)
# Output:   - dip_ksi.txt -> 1 to 1 relation between a kinase and a substrate. In database format: source, target, source_name...with three more columns: evidence (evidence codes), interaction_id (in the file), expanded (if more than one kinase was present)
#           - dip_ppi.txt 
###################################################################################################################################
import argparse
import re
import copy
from xml.etree import ElementTree as et
import mappings as map
import datetime


today = datetime.date.today().strftime("%Y%m%d")
# parse arguments
parser = argparse.ArgumentParser(description="Extraction of KSI from PID XML files.")

parser.add_argument("pidfile")
parser.add_argument("-e", "--expanded", action="store_true", help="consider KSI when several kinases involved (one KSI per kinase)")
parser.add_argument("-m", "--members", action="store_true", help="consider family memebers")
parser.add_argument("-o", "--outputfile", help="name of output file")

args = parser.parse_args()

# create list of kinases uniprot accession codes
kinacc2name = map.getKinaseAcc2Symbol()

# get dictionary of entrez to uniprot mainly for BioCarta
entrez2uniprot = map.getEntrez2Uniprot()

# get dictionary of uniprot accession codes to gene symbols
uniprot2symbol = map.getUniprotMapDicts()
sec2prim = map.getSecondary2PrimaryAcc()

# read PID file
print("\nReading PID XML file...\n")

up_regex = "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"

tree=et.parse(open(args.pidfile))
root=tree.getroot()

molecules = root.find("Model").find("MoleculeList").findall("Molecule")# list of molecules

id2uniprot = {}
entrez_nomap = []
print("\tdoing PROTEINS...")

# first extract simple proteins (not complexes) -> necessary to then map complexes and proteins with "members"
for molecule in molecules:
    mol_type = molecule.attrib["molecule_type"]
    
    if mol_type == "protein": # should be only one component
        
        mol_id = molecule.attrib["id"]
        #print("\t\tdoing protein of id %s"%(mol_id))
        mol_up = []
        mol_pf = []
        
        for name in molecule.iter("Name"):
            name_type = name.attrib["name_type"]
            if name_type == "UP": # uniprot
                name_value = name.attrib["value"].split("-")[0] #isoform
                if re.match(up_regex, name_value):
                    mol_up.append(name_value) 
            elif name_type == "ZZ": # no idea but sometimes this also contains uniprot codes...but with isoform notation i.e. ".2" etc
                name_value = name.attrib["value"].split(".")[0].upper()
                name_value = name_value.split("-")[0] #isoform
                if re.match(up_regex, name_value):
                    mol_up.append(name_value)
            elif name_type == "LL":
                name_value = name.attrib["value"] # entrez gene (BioCarta)
                try:
                    name_value = entrez2uniprot[name_value]
                    mol_up.append(name_value)
                except:
                    entrez_nomap.append(name_value)
            elif name_type == "PF" or "OF": # preferred or official symbol (sometimes both exist...)
                name_value = name.attrib["value"] # preferred symbol
                mol_pf.append(name_value)

        if mol_up != []:
            id2uniprot[mol_id] = {}
            
            id2uniprot[mol_id]["name"] = mol_pf 
            id2uniprot[mol_id]["uniprot"] = mol_up # can be uniprot or entrez (for BioCarta, mainly)
            id2uniprot[mol_id]["ptms"] = [[]] # individual proteins don't have ptm

ppi_dic = {}
ppi_dic = {}

# apply filter for members
if args.members:
    moltypes = ['protein', 'complex']
else:
    moltypes = ['complex']

for moltype in moltypes: #first parse the "protein" Families, then the "complexes" (only possible option, as some complexes have families)
    print("\tdoing %s..."%(re.sub("PROTEIN", "PROTEIN FAMILIES", moltype.upper())))
    if moltype == "complex":
        iter_name = "ComplexComponent"
        iter_att = "molecule_idref"
    elif moltype == "protein":
        iter_name = "Member"
        iter_att = "member_molecule_idref"
   
    for molecule in molecules:
        complex_elements = []
        mol_type = molecule.attrib["molecule_type"]

        if mol_type == moltype:
            
            if mol_type == "protein" and molecule.find("FamilyMemberList") is None: # not a Family
                continue

            mol_id = molecule.attrib["id"]
            mol_name = molecule.find("Name").attrib["value"] 
            
            #print("\t\tdoing %s of id %s"%(moltype, mol_id))
            
            mol_up_ptm = {}
            for component in molecule.iter(iter_name):
                comp_id = component.attrib[iter_att]
                if comp_id in id2uniprot.keys():
                    comp_up = copy.deepcopy(id2uniprot[comp_id]["uniprot"])
                    
                    # append complex elements
                    if mol_type == "complex":
                        complex_elements.extend(comp_up)
                    
                    # annotate ptms
                    for comp in comp_up:
                        mol_up_ptm[comp] = []
                        
                    for ptm in component.iter("PTMTerm"):
                        ptm_protein = ptm.attrib["protein"].split("-")[0]
                        if ptm_protein == "":
                            if len(comp_up) == 1:
                                ptm_protein = comp_up[0]
                            else:
                                continue
                        ptm_position = ptm.attrib["position"]
                        ptm_mod = ptm.attrib["modification"]

                        ptm_whole = ptm_position + "||" + ptm_mod 
                        if ptm_whole:
                            mol_up_ptm[ptm_protein].append(ptm_whole)
                        else:
                            mol_up_ptm[ptm_protein].append("")
        
            if mol_up_ptm != {}:
                id2uniprot[mol_id] = {}
            
                id2uniprot[mol_id]["name"] = mol_name
                id2uniprot[mol_id]["uniprot"] = mol_up_ptm.keys()
                id2uniprot[mol_id]["ptms"] = mol_up_ptm.values()
            
        if complex_elements != []:
            for i in range(0, len(complex_elements)-1):
                for j in range(i+1, len(complex_elements)):
                    pair = [complex_elements[i], complex_elements[j]]
                    pair.sort()
                    key = '||'.join(pair)

                    if key in ppi_dic.keys():
                        if mol_id not in ppi_dic[key]["complex_id"]:
                            ppi_dic[key]["complex_id"].append(mol_id)
                            ppi_dic[key]["complex_nb"].append(str(len(complex_elements)))
                            ppi_dic[key]["complex_name"].append(mol_name)
                    else: 
                        ppi_dic[key] = {}

                        ppi_dic[key]["source"] = pair[0]
                        ppi_dic[key]["target"] = pair[1]
                        ppi_dic[key]["complex_id"] = [mol_id]
                        ppi_dic[key]["complex_nb"] = [str(len(complex_elements))]
                        ppi_dic[key]["complex_name"] = [mol_name]
                    

# retrieve phosphorylations (only for PID or Reactome, not BioCarta - ambiguous) 
if not re.search("BioCarta", args.pidfile):
    interactions = root.find("Model").find("InteractionList")
    print("\tdoing INTERACTIONS...")
else:
    interactions = [] # empty list in case of interactions

ksi_dic = {}
for interaction in interactions:
    int_type = interaction.attrib["interaction_type"]
    int_id_all = interaction.attrib["id"]
     
    # phosphorylation witin type "modification"
    if int_type == "modification": 
        
        #print("\t\tdoing interaction with id %s"%(int_id_all))
        # get inputs, outputs, and agents (catalyzer)
        inputs = {}
        outputs = {}
        agents = []
        for interactor in interaction.iter("InteractionComponent"):
            int_role = interactor.attrib["role_type"]
            int_id = interactor.attrib["molecule_idref"]
        
            try:
                ints = dict(zip(copy.deepcopy(id2uniprot[int_id]["uniprot"]),
                                copy.deepcopy(id2uniprot[int_id]["ptms"]))) 
            except:
                continue
                    
            for ptm in interactor.iter("PTMTerm"):
                ptm_protein = ptm.attrib["protein"].split("-")[0]
                ptm_position = ptm.attrib["position"]
                ptm_mod = ptm.attrib["modification"]
                
                ptm_whole = ptm_position + "||" + ptm_mod
                    
                if ptm_protein in ints.keys() and ptm_whole not in ints[ptm_protein]:
                    ints[ptm_protein].append(ptm_whole)
                else:
                    print("\tReaction %s: ptm protein %s not in %s"%(int_id_all, ptm_protein, ','.join(ints.keys())))
            
            if int_role == "input":
                for int in ints.keys():
                    if int not in inputs.keys():
                        inputs[int] = ints[int]
                    else:
                        if ints[int] not in inputs[int]:
                            inputs[int].extend(ints[int]) 
            elif int_role == "output":
                for int in ints.keys():
                    if int not in outputs.keys():
                        outputs[int] = ints[int]
                    else:
                        if ints[int] not in outputs[int]:
                            outputs[int].extend(ints[int])
            elif int_role == "agent":
                agents.extend(ints.keys())

        # phosphorylation needs the three elements: input, output, agent
        if not inputs or not outputs or not agents:
            continue
    
        # check that agent(s) is kinase -> if several proteins, take the kinase
        kinase = []
        for agent in agents:        
            if agent in kinacc2name.keys():
                kinase.append(agent)
        kinase = list(set(kinase))
        
        if kinase == []:
            continue

        # apply filter for kinase "expansion"
        if len(kinase)>1 and not args.expanded:
            continue

        # check phosphorylation event: compare inputs vs. outputs
        for input in inputs.keys():
            prot = input
            input_mod = inputs[prot]
            try:
                output_mod = outputs[prot]
                #out_index = outputs.index(prot) # check input UniProt is in output
            except:
                print("\tReaction %s: input %s not in output %s"%(int_id_all, prot, ','.join(outputs)))
                continue
            #output_mod = outputs_mods[out_index]
            mod_diff = list(set(output_mod).difference(input_mod))
            phospho_positions = []
            for mod in mod_diff:
                if re.search("phosphorylation", mod):
                    phosph = mod.split("||")[0]
                    phospho_positions.append(phosph)
            
            if phospho_positions != []:
                #print("\nKSI (%s):"%(int_id_all))
                #print("\tkinase: %s"%','.join(kinase))
                #print("\tsubstrate %s"%(prot))
                #print("\tmodification %s"%(','.join(mod_diff)))
                
                pmids = []
                for ref in interaction.iter("Reference"):
                    pmid = ref.attrib["pmid"]
                    pmids.append(pmid)

                evids = []
                for ev in interaction.iter("Evidence"):
                    evid = ev.attrib["value"]
                    evids.append(evid)
               
                # if more than one kinase, duplicate the entry (should it be like this???)
                if len(kinase) >1:
                    kin_exp = 'yes'
                else:
                    kin_exp = 'no'

                for kin in kinase:
                    key = kin + "||" + prot

                    if not ksi_dic.has_key(key):
                        ksi_dic[key] = {}

                        ksi_dic[key]['source'] = kin
                        ksi_dic[key]['target'] = prot
                        ksi_dic[key]['pmids'] = pmids
                        ksi_dic[key]['evids'] = evids
                        ksi_dic[key]['positions'] = phospho_positions
                        ksi_dic[key]['id'] = int_id_all
                        ksi_dic[key]['expanded_kinases'] = kin_exp
                    else:
                        
                        for pmid in pmids:
                            if pmid not in ksi_dic[key]['pmids']:
                                ksi_dic[key]['pmids'].append(pmid)
                        for evid in evids:
                            if evid not in ksi_dic[key]['evids']:
                                ksi_dic[key]['evids'].append(evid)
                        for phosph in phospho_positions:
                            if phosph not in ksi_dic[key]['positions']:
                                ksi_dic[key]['positions'].append(phosph)
            
# write files:
nomap=[]            
print("\n Writing files...\n")

if args.outputfile:
    outfilename = args.outputfile
else:
    outfilename = "pid_" + args.pidfile.split(".")[0]
    if args.members:
        outfilename += "_members"
    if args.expanded:
        outfilename += "_expanded"
    outfilename += "_ksi" + today + ".txt"

outfields = ["source", "target", "source_name", "target_name", "PMIDs", "dates", "sources", "type", "positions"]
outfields.extend(["evidence", "interaction_id"])
if args.expanded:
    outfields.append("expanded_kinases")

dates = "2014"
sources = "pid"#(" + args.pidfile.split(".")[0] + ")"
type = "KSI"
with open(outfilename, "w") as outfile:

    for field in outfields:
        outfile.write(field+"\t")
    outfile.write("\n")

    for key in ksi_dic.keys():
        source = ksi_dic[key]['source']
        target = ksi_dic[key]['target']
            
        try:
            source_name = uniprot2symbol[source]
        except:
            try:
                source = sec2prim[source]
                source_name = uniprot2symbol[source]
            except:
                nomap.append(source)
                continue

        try:
            target_name = uniprot2symbol[target]
        except:
            try:
                target = sec2prim[target]
                target_name = uniprot2symbol[target]
            except:
                nomap.append(target)
                continue

        PMIDs = ','.join(ksi_dic[key]['pmids'])

        interaction_id = ksi_dic[key]['id']
            
        positions = ",".join(ksi_dic[key]['positions'])
            
        evidence = ",".join(ksi_dic[key]['evids'])
        
        expanded_kinases = ksi_dic[key]['expanded_kinases']
        for field in outfields:
            outfile.write(vars()[field] + "\t")
        outfile.write("\n")

print("\nNot mapped accessions for KSI: %d being %s"%(len(set(nomap)), ','.join(set(nomap))))

# write second file:
outfields = ["source", "target", "source_name", "target_name", "PMIDs", "dates", "sources", "type", "source_is_bait", "target_is_bait"]
outfields.extend(["complex_name", "complex_id", "complex_nb"])
nomap_ppi = []
PMIDs = ""
dates = "2014"
sources = "pid"#(" + args.pidfile.split(".")[0] + ")"
type = "PPI"
source_is_bait = "no"
target_is_bait = "no"

outfilename = "pid_" + args.pidfile.split(".")[0]
if args.members:
    outfilename += "_members"
outfilename += "_ppi_" + today + ".txt"

with open(outfilename, "w") as outfile:
    for field in outfields:
        outfile.write(field + "\t")
    outfile.write("\n")
    
    for key in ppi_dic.keys():
        source = ppi_dic[key]["source"]
        target = ppi_dic[key]["target"]

        try:
            source_name = uniprot2symbol[source]
        except:
            try:
                source = sec2prim[source]
                source_name = uniprot2symbol[source]
            except:
                nomap_ppi.append(source)
                continue
        try:
            target_name = uniprot2symbol[target]
        except:
            try:
                target = sec2prim[target]
                target_name = uniprot2symbol[target]
            except:
                nomap_ppi.append(target)
                continue
        
        complex_name = '~'.join(ppi_dic[key]["complex_name"])
        complex_id = ','.join(ppi_dic[key]["complex_id"])
        complex_nb = ','.join(ppi_dic[key]["complex_nb"])

        for field in outfields:
            outfile.write(vars()[field] + "\t")
        outfile.write("\n")

print("\nNot mapped accessions for PPI: %d being %s"%(len(set(nomap_ppi)), ','.join(set(nomap_ppi))))
print("\nNot mapped Entrez Gene ids to UniProt: %d being %s"%(len(set(entrez_nomap)),','.join(set(entrez_nomap))))
