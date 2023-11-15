#!/usr/bin/env python
# Adapted from souporcell_pipeline.py

import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--bam", required=True, type=str, help="")
parser.add_argument("--n_splits", required=True, type=int, help="")
parser.add_argument("--out", required=True, type=str, help="The file where results will be saved.")
args = parser.parse_args()

print(os.path.dirname(args.out))
if not os.path.isdir(os.path.dirname(args.out)):
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

import pysam
import math
import os
import json


print("Reading BAM file.")
bam = pysam.AlignmentFile(args.bam)
total_reference_length = 0
for chrom in bam.references:
    total_reference_length += bam.get_reference_length(chrom)
step_length = int(math.ceil(total_reference_length / args.n_splits))
print("\tStep length: {:,}.".format(step_length))

regions = {}
region = []
region_so_far = 0
chrom_so_far = 0
index = 0
for chrom in bam.references:
    chrom_length = bam.get_reference_length(chrom)
    while True:
        if (chrom_length - chrom_so_far) <= step_length - region_so_far:
            region.append({"chr": chrom, "start": chrom_so_far, "stop": chrom_length})
            region_so_far += chrom_length - chrom_so_far
            chrom_so_far = 0
            break
        else:
            region.append({"chr": chrom, "start": chrom_so_far, "stop": chrom_so_far + step_length - region_so_far})
            regions[index] = region
            index += 1
            region = []
            chrom_so_far += step_length - region_so_far
            region_so_far = 0

if len(region) > 0:
    regions[index] = region

print("Regions:")
for index, region_list in regions.items():
    print("\t{}: {} regions".format(index, len(region_list)))
print("")

if len(regions) != args.n_splits:
    print("Error, created more splits than exepcted.")
    exit()

print("Saving {}".format(os.path.basename(args.out)))
with open(args.out, 'w') as f:
    json.dump(regions, f, indent=4)
f.close()

print("Done.")
