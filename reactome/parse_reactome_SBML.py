#!/cm/shared/apps/python/2.7.6/bin/python

###################################################################################################################################
# Script to retrieve KSI from Reactome SBML file "Human reactions in SBML" at reactome.org/download/
# 
# Uses as input: - reactome SBML file
#
# Checks if it's phosphorylation by looking at ATP as reactant and ADP as product.
# Filters kinases so that we only include the ones in the kinase list 
#
# Output:   - reactome_ksi.txt -> 1 to 1 relation between a kinase and a substrate. In database format: source, target, source_name...
#           - reactome_ksi_ambiguous -> many to many relation between kinases and substrates. Requires manual curation to know which
#                                       kinase phosphorylates which substrate. Format like in database with added columns: reaction_id, 
#                                       reaction_name, description, etc...to help curation.
###################################################################################################################################
import argparse
import re
from xml.etree import ElementTree as et
import mappings as map

# parse arguments
parser = argparse.ArgumentParser(description="Extraction of KSI from Reactome SBML files.")

parser.add_argument("reactomefile")

args = parser.parse_args()

# create list of kinases uniprot accession codes
uniprot2symbol, syn2uniprot = map.getUniprotMapDicts()

# create dictionary of uniprot accession codes to gene symbols (if no symbol:empty string)
kinacc2name = map.getKinaseAcc2Symbol()

# read reactome file
print("\nReading Reactome SBML file...\n")

up_regex = "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"

tree = et.parse(open(args.reactomefile))
root = tree.getroot()

# map species names and uniprot codes
spid2uniprot = {}
spid2name = {}
key2pmid = {}
for sp in root.getiterator("{http://www.sbml.org/sbml/level2/version4}species"):
    uniprots = []
    
    sp_name = sp.attrib["name"]
    sp_id = sp.attrib["id"]
    sp_metaid = sp.attrib["metaid"]
    sp_sbo = sp.attrib["sboTerm"]
    #print('\n%s: '%(sp_name))

    for itis in sp.getiterator("{http://biomodels.net/biology-qualifiers/}is"):
        #print('\tis:') #it should only contain 1...
        for li in itis.getiterator("{http://www.w3.org/1999/02/22-rdf-syntax-ns#}li"):
            sp_acc = li.attrib['{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource'].split(":")[-1]
            sp_acc = sp_acc.split("-")[0] #in case of isoform naming i.e. P12345-2
            if re.match(up_regex, sp_acc):
                uniprots.append(sp_acc)
                #print('\t\t%s'%(sp_acc))
        
    for haspart in sp.getiterator("{http://biomodels.net/biology-qualifiers/}hasPart"):
        #print('\thas part\n')
        for li in haspart.getiterator("{http://www.w3.org/1999/02/22-rdf-syntax-ns#}li"):
            sp_acc = li.attrib['{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource'].split(":")[-1]
            sp_acc = sp_acc.split("-")[0]
            if re.match(up_regex, sp_acc):
                uniprots.append(sp_acc)
                #print('\t\t%s'%(sp_acc))
    
    spid2name[sp_id] = sp_name
    if len(uniprots) >= 0:
        spid2uniprot[sp_id] = uniprots

#print(spid2uniprot)

print("\nExtracting KSI...\n")
ksi_dic = {}
kin_out = [] # kinases not found in my list from uniprot
for react in root.getiterator("{http://www.sbml.org/sbml/level2/version4}reaction"):
    react_name = react.attrib["name"]
    react_id = react.attrib["id"]
    for desc in react.getiterator("{http://www.w3.org/1999/xhtml}p"):
        react_desc = desc.text
    #if not re.search(" phosphorylates ", react_name):
    #    continue
    #print("\tREACTION: %s"%(react_name))
    for reactant in react.getiterator("{http://www.sbml.org/sbml/level2/version4}listOfReactants"):
        flag = 0 # to check that it's a real phosphorylation   
        subs = []
        subs_names = []
        for sp in reactant.getiterator("{http://www.sbml.org/sbml/level2/version4}speciesReference"):
            sp_id = sp.attrib["species"]
            
            if re.match("ATP", spid2name[sp_id].split(" ")[0]):
                flag = 1
            else:
                sp_acc = spid2uniprot[sp_id]
                    
                if sp_acc != []:
                    subs.extend(sp_acc)
                
                sp_name = spid2name[sp_id]
                subs_names.append(sp_name)
    
    if flag == 0 or subs == []:
        continue

#    print("\tREACTION: %s"%(react_name))
#    print("\t\tSubstrate: %s"%(','.join(subs)))
    for catalyser in react.getiterator("{http://www.sbml.org/sbml/level2/version4}listOfModifiers"):
        kin = []
        kin_names = []
        for sp in catalyser.getiterator("{http://www.sbml.org/sbml/level2/version4}modifierSpeciesReference"):
            sp_id = sp.attrib["species"]

            sp_acc = spid2uniprot[sp_id]

            if sp_acc != []:
                kin.extend(sp_acc)
            
            sp_name = spid2name[sp_id]
            kin_names.append(sp_name)

    # if more than one kinase (because of a complex for example), check which of the components of the complex is the actual kinase and keep that one
    kin_filt = []
    for k in kin:
        if k in kinacc2name.keys():
            kin_filt.append(k)
        else:
            kin_out.append(k)
                    
    if kin_filt == []:
        continue

    for prod in react.getiterator("{http://www.sbml.org/sbml/level2/version4}listOfProducts"):
        prods = []
        prods_names = []
        flag = 0
        for sp in prod.getiterator("{http://www.sbml.org/sbml/level2/version4}speciesReference"):
            sp_id = sp.attrib['species']
            
            if re.match("ADP", spid2name[sp_id].split(" ")[0]):
                flag = 1
            else:
                sp_acc = spid2uniprot[sp_id]
                if sp_acc != []:
                    prods.extend(sp_acc)
                
                sp_name = spid2name[sp_id]
                prods_names.append(sp_name)

    if flag == 0:
        continue
#    print("\t\tKinase: %s"%(','.join(kin)))
    # get PMIDs
    pmids = []
    for desc in react.getiterator("{http://biomodels.net/biology-qualifiers/}isDescribedBy"):
        for li in desc.getiterator("{http://www.w3.org/1999/02/22-rdf-syntax-ns#}li"):
            pmid = li.attrib['{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource'].split(":")[-1]
            pmids.append(pmid)
    
    # add to dictionary
    key = ','.join(kin) + '||' + ','.join(subs)
    if key not in ksi_dic.keys():
        ksi_dic[key] = {}

        ksi_dic[key]['source'] = kin_filt
        ksi_dic[key]['target'] = subs
        ksi_dic[key]['product'] = prods
        ksi_dic[key]['pmids'] = pmids
        ksi_dic[key]['react_name'] = react_name
        ksi_dic[key]['react_desc'] = react_desc
        ksi_dic[key]['react_id'] = react_id
        ksi_dic[key]['src_name'] = kin_names
        ksi_dic[key]['tar_name'] = subs_names
        ksi_dic[key]['prod_name'] = prods_names
        
    else:
        ksi_dic[key]['pmids'].extend(pmids)

# write file
print("\nWriting files...\n")
outfilename1 = "reactome_ksi.txt"
outfilename2 = "reactome_ksi_ambiguous.txt"
outfields = ["source", "target", "source_name", "target_name", "PMIDs", "dates", "sources", "type", "phospho_positions"]
    
dates = "2014"
sources = "reactome"
type = "KSI"
phospho_positions = ""
nomap = []
with open(outfilename1, "w") as outfile1, open(outfilename2, "w") as outfile2:
        for field in outfields:
            outfile1.write(field+"\t")
            outfile2.write(field+"\t")
        outfile2.write("reaction_name" + "\t" + "reaction_id" + "\t" + "source_orig_names" + "\t" + "tar_orig_names" + "\t" + "prods_orig_names" + "\t" + 
                       "product" + "\t" + "product_name" + "\t" + "full_description" + "\t")
        outfile1.write("\n")
        outfile2.write("\n")

        for key in ksi_dic.keys():
            source_all = ksi_dic[key]['source']
            target_all = ksi_dic[key]['target']
            product_all = ksi_dic[key]['product']
            
            source_name = []
            target_name = []
            product_name = []
            source = []
            target = []
            product = []
            for src in source_all:
                if src in uniprot2symbol.keys():
                    src_name = uniprot2symbol[src]
                    source_name.append(src_name)
                    source.append(src)
                else:
                    nomap.append(src)
            for tar in target_all:
                if tar in uniprot2symbol.keys():
                    tar_name = uniprot2symbol[tar]
                    target_name.append(tar_name)
                    target.append(tar)
                else:
                    nomap.append(tar)
            for prod in product_all:
                if prod in uniprot2symbol.keys():
                    prod_name = uniprot2symbol[prod]
                    product_name.append(prod_name)
                    product.append(prod)
            
            # choose file
            if len(source) == 1 and len(target) == 1:
                outfile = outfile1
            else:
                outfile = outfile2

            source = ','.join(source)
            target = ','.join(target)
            source_name = ','.join(source_name)
            target_name = ','.join(target_name)
            product_name = ','.join(product_name)

            PMIDs = ','.join(ksi_dic[key]['pmids'])

            react_name = ksi_dic[key]['react_name']
            react_id = ksi_dic[key]['react_id']
            react_desc = ksi_dic[key]['react_desc'].encode("utf-8",'ignore')
            prods = ",".join(ksi_dic[key]['product'])
            
            source_names = ",".join(ksi_dic[key]['src_name'])
            tar_names = ",".join(ksi_dic[key]['tar_name'])
            prods_names = ",".join(ksi_dic[key]['prod_name'])
            
            
            for field in outfields:
                outfile.write(vars()[field] + "\t")
                
            if outfile == outfile2:
                outfile.write(react_name + "\t" + react_id + "\t" + source_names + "\t" + tar_names + "\t" + 
                              prods_names + "\t" + prods + "\t" + product_name + "\t" + react_desc + "\t")
            outfile.write("\n")

print("\nKinases not in kinase list: %d being %s"%(len(set(kin_out)), ','.join(set(kin_out))))

print("\nProteins not mapped to UniProt: %d being %s"%(len(set(nomap)), ','.join(set(nomap))))
