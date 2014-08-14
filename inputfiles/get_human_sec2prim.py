#!/cm/shared/apps/python/2.7.6/bin/python
######################################################################################
# Select only the human secondary to primary identifiers from the whole uniprot list.
# Splitting cases are not considered (where a secondary accession gives two primary)
######################################################################################
import mappings as map
import argparse

# parse arguments
parser = argparse.ArgumentParser(description="Selection of human secondary to primary uniprot accessions")
parser.add_argument("sec2primfile")


# get list of human uniprot accessions
print("\nReading list of human uniprot primar accessions...\n")
up2name = map.getUniprotMapDicts()
args = parser.parse_args()

infilename = args.sec2primfile

# select only human sec2prim
sec2prim = {}
reps = []
with open(infilename) as infile:
    print("\nSelecting human uniprot secondary accession codes to primary accession codes...\n")

    lines = infile.readlines()
    indx=0
    for line in lines:
        indx += 1
        if indx%10000 == 0:
            print('\t%d'%indx)


        split_line = line.rstrip().split("\t")

        sec_acc = split_line[0]
        prim_acc = split_line[1]

        if prim_acc in up2name.keys(): #if its human
            if sec_acc in sec2prim.keys():
                reps.append(sec_acc)
            else:
                sec2prim[sec_acc] = prim_acc
    reps = set(reps)
    for rep in reps:
        del sec2prim[rep]

# print file
print("\nWriting file...\n")
outfilename = infilename.split(".")[0] + "_human.txt"
with open(outfilename, 'w') as outfile:
    for key in sec2prim.keys():
        outfile.write(key + "\t" + sec2prim[key] + "\n")
