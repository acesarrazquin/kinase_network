#!/cm/shared/apps/python/2.7.6/bin/python

##########################################
#   Reorder columns
#
##########################################

import sys


headers = ["source", "target", "source_name", "target_name", "PMIDs", "dates", "sources", "type"]
headers.extend(['source_IDs', 'target_form', 'strength', 'strength_type',
                 'strength_direction', 'strong_target', 'technology', 'cell_type'])

filename = sys.argv[1]
outfilename = filename + 'corrected'
with open(filename) as infile, open(outfilename, 'w') as outfile:
    lines = infile.readlines()

    header = lines[0]
    fields = header.rstrip().split('\t')
    for header in headers:
        outfile.write(header + '\t')
    outfile.write("\n")

    for line in lines[1:len(lines)]:

        split_line = line.split("\n")[0].split('\t')


        headers = ["source", "target", "source_name", "target_name", "PMIDs", "dates", "sources", "type"]
        headers.extend(['source_IDs', 'target_form', 'strong_target', 'strength', 'strength_type',
                'strength_direction', 'technology', 'cell_type'])

        for i in range(0, len(headers)):
            vars()[headers[i]] = split_line[i]

        type="DPI"
        headers = ["source", "target", "source_name", "target_name", "PMIDs", "dates", "sources", "type"]
        headers.extend(['source_IDs', 'target_form', 'strength', 'strength_type',
                        'strength_direction', 'strong_target', 'technology', 'cell_type'])
        for header in headers:
            print(header)
            value = vars()[header]
            if value==" ":
                value=""

            outfile.write(value + "\t")
        outfile.write("\n")


