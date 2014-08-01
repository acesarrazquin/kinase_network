#!/cm/shared/apps/python/2.7.6/bin/python

################################################################################################################
# compile_interactions.py takes PPI/KSI/DPI files and merges identical source-target interactions, keeping track
# of sources, years, etc...
#
# If KSI and PPI are equal, then it's considered as KSI.
################################################################################################################
import argparse

parser = argparse.ArgumentParser(description="Assembly of KSI from NetworKIN predictions.")

parser.add_argument("files", nargs="+")

args = parser.parse_args()

print(args)