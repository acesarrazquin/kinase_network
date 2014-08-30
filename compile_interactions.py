#!/cm/shared/apps/python/2.7.6/bin/python

################################################################################################################
# compile_interactions.py takes PPI/KSI/DPI files and merges identical source-target interactions, keeping track
# of sources, years, etc...
#
# Note: If KSI and PPI are equal, then it's considered as KSI.
#
# Options: select type of interaction, remove self-interactions...
################################################################################################################
import argparse
import re
import datetime


def appendKey(int_dic, key, headers, values):
    int_dic[key] = {}
    for i in range(0, len(headers)):
        key2 = headers[i]
        value = values[i]
        int_dic[key][key2] = value

def updateKey(int_dic, key, **kw):
    for key2 in kw.keys():
        if kw[key2] == '' or kw[key2] == "NA":
            continue

        if int_dic[key][key2] != "" and int_dic[key][key2] != "NA":
            int_dic[key][key2] += ',' + kw[key2]
        else:
            int_dic[key][key2] = kw[key2]

################### MAIN####################################

parser = argparse.ArgumentParser(description="Assembly of KSI from NetworKIN predictions.")

parser.add_argument("files", nargs="+") # variable number of arguments (file names)
parser.add_argument("--type", "-t", choices=["PPI", "KSI", "DPI", "all"], dest="int_type", default="all", help="Type of interaction")
parser.add_argument("--no_self", "-s", action="store_true", help="Don't consider self-interactions", dest="noself")
parser.add_argument("--output_name", "-o", dest="outname", help="Name of output file")
parser.add_argument("--dpi_strong", "-d", action ="store_true", dest="strong", help="only strong DPI")

args = parser.parse_args()

# choose output headers and filenames
headers = ["source", "target", "source_name", "target_name", "PMIDs", "dates", "sources", "type"]
today = datetime.date.today().strftime("%Y%m%d")
if args.int_type == "PPI":
    headers.extend(["source_is_bait", "target_is_bait"])
    outname = "all_ppi_%s.txt"%(today)
elif args.int_type == "KSI":
    headers.append("positions")
    outname = "all_ksi_%s.txt"%(today)
elif args.int_type == "DPI":
    headers.extend(['source_IDs', 'target_form', 'strength', 'strength_type',
                    'strength_direction', 'strong_target', 'technology', 'cell_type'])
    if args.strong:
        outname = "all_strong_dpi_%s.txt"%(today)
    else:
        outname = "all_dpi_%s.txt"%(today)
else:
    outname = "full_network_%s.txt"%(today)

# Read files
print("\nReading interaction files...\n")

as_ksi = []
int_dic = {}
for filename in args.files:
    with open(filename) as infile:
        lines = infile.readlines()

        for line in lines[1:len(lines)]:

            # get variables
            values = line.split("\n")[0].split("\t")
            for i in range(0, len(headers)):
                vars()[headers[i]] = values[i]

            # skip if type doesn't match. Note: This is problematic because "type" is a reserved word in Python. Now it's overwritten.
            if args.int_type != "all" and type != args.int_type:
                continue

            # for DPI: if only strong required, skip if not strong
            if args.int_type == "DPI" and args.strong and strong_target !="yes":
                continue

            # apend to dictionary
            key = "%s|%s"%(values[0], values[1])
            if not int_dic.has_key(key):
                appendKey(int_dic, key, headers, values)

            else:
                if int_dic[key]["type"] == type:
                    if args.int_type == "PPI":
                        updateKey(int_dic, key, PMIDs=PMIDs, sources=sources, dates=dates,
                                  source_is_bait=source_is_bait, target_is_bait=target_is_bait)
                        for bait in ["source_is_bait", "target_is_bait"]:
                            if re.search("yes", int_dic[key][bait]):
                                int_dic[key][bait] = "yes"
                            else:
                                int_dic[key][bait] = "no"
                    elif args.int_type == "KSI":
                        updateKey(int_dic, key, PMIDs=PMIDs, sources=sources, dates=dates,
                                  positions=positions)
                    elif args.int_type == "DPI":
                        updateKey(int_dic, key, PMIDs=PMIDs, sources=sources, dates=dates,
                                  target_form=target_form, strength=strength, strength_type=strength_type,
                                  strength_direction=strength_direction, strong_target=strong_target,
                                  technology=technology, cell_type=cell_type)
                    else:
                        updateKey(int_dic, key, PMIDs=PMIDs, sources=sources, dates=dates)

                else:
                    # consider only the KSI; if new is KSI overwrite the previous
                    if int_dic[key]["type"] == "PPI" and type == "KSI": # it's redundant normally
                        #print("\t%s will be only considered as KSI..."%(key))
                        as_ksi.append(key)
                        appendKey(int_dic, key, headers, values)
                    else:
                        print("\t%s error! %s and %??!!"%(key, int_dic[key]["type"], type))

# write the final file
print("\nWriting final file...\n")

if args.outname:
    outfilename = args.outname
else:
    outfilename = outname

headers_set = ["PMIDs", "sources",
               "positions",
               "source_IDs"]

with open(outfilename, "w") as outfile:

    outfile.write('\t'.join(headers) + '\n')

    for key in int_dic.keys():
        values = ""
        for header in headers:
            if header == "dates":
                value = ",".join(sorted(set(int_dic[key][header].split(",")), reverse=True))
            elif header in headers_set:
                value = ",".join(sorted(set(int_dic[key][header].split(","))))
            else:
                value = int_dic[key][header]
            values += "\t" + value
        values = values[1:len(values)] + "\n"
        outfile.write(values)

as_ksi = set(as_ksi)
print("\nSUMMARY: considered as KSI %d, being %s" %(len(as_ksi), ','.join(as_ksi)))
