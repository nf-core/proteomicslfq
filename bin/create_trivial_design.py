#!/usr/bin/env python3

import sys
import glob
import re

# code to sort in a human readable way
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    e.g. UPS1_50amol_R1 is less than UPS1_1200amol_R1
    http://nedbatchelder.com/blog/200712/human_sorting.html
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

if not len(sys.argv) == 3:
  print("Usage: MZML_FOLDER LABEL_PER_FILE")
  exit()

in_path = sys.argv[1]
label_per_file = int(sys.argv[2])

mzmls = [f for f in glob.glob(in_path + "/*.mzML", recursive=False)]
mzmls.sort(key=natural_keys)

file_count = 1
fraction_group = 1
label = 1
sample = 1
print("Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample")

for f in mzmls:
  for label in range(1, label_per_file + 1): 
    print(str(file_count) + "\t" + str(fraction_group) + "\t" + f + "\t" + str(label) + "\t" + str(sample))
    sample += 1
  file_count += 1
print()
print("Sample\tMSstats_Condition\tMSstats_BioReplicate")
sample = 1
condition = 1
bioreplicate = 1
for f in mzmls:
  for label in range(1, label_per_file + 1): 
    print(str(sample) + "\t" + str(condition) + "\t" + str(bioreplicate))
    sample += 1
    condition += 1
    bioreplicate += 1
