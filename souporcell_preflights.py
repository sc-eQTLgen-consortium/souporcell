#!/usr/bin/env python
# Adapted from souporcell_pipeline.py

import argparse
import os
import json

parser = argparse.ArgumentParser(
    description="single cell RNAseq mixed genotype clustering using sparse mixture model clustering.")
parser.add_argument("-i", "--bam", required = True, help = "cellranger bam")
parser.add_argument("-b", "--barcodes", required = True, help = "barcodes.tsv from cellranger")
parser.add_argument("-f", "--fasta", required = True, help = "reference fasta file")
parser.add_argument("--n_remap_splits", required = True, type = int, help = "")
parser.add_argument("--n_freebayes_splits", required = True, type = int, help = "")
parser.add_argument("-o", "--out_dir", required = True, help = "name of directory to place souporcell files")
parser.add_argument("-k", "--clusters", required = True, help = "number cluster, tbd add easy way to run on a range of k")
parser.add_argument("-p", "--ploidy", required = False, default = "2", help = "ploidy, must be 1 or 2, default = 2")
parser.add_argument("--min_alt", required = False, default = "10", help = "min alt to use locus, default = 10.")
parser.add_argument("--min_ref", required = False, default = "10", help = "min ref to use locus, default = 10.")
parser.add_argument("--max_loci", required = False, default = "2048", help = "max loci per cell, affects speed, default = 2048.")
parser.add_argument("--restarts", required = False, default = 100, type = int,
    help = "number of restarts in clustering, when there are > 12 clusters we recommend increasing this to avoid local minima")
parser.add_argument("--common_variants", required = False, default = None,
    help = "common variant loci or known variant loci vcf, must be vs same reference fasta")
parser.add_argument("--known_genotypes", required = False, default = None,
    help = "known variants per clone in population vcf mode, must be .vcf right now we dont accept gzip or bcf sorry")
parser.add_argument("--known_genotypes_sample_names", required = False, nargs = '+', default = None,
    help = "which samples in population vcf from known genotypes option represent the donors in your sample")
parser.add_argument("--skip_remap", required = False, default = False, type = bool,
    help = "don't remap with minimap2 (not recommended unless in conjunction with --common_variants")
parser.add_argument("--no_umi", required = False, default = "False", help = "set to True if your bam has no UMI tag, will ignore/override --umi_tag")
parser.add_argument("--umi_tag", required = False, default = "UB", help = "set if your umi tag is not UB")
parser.add_argument("--cell_tag", required = False, default = "CB", help = "DOES NOT WORK, vartrix doesnt support this! set if your cell barcode tag is not CB")
args = parser.parse_args()

if not os.path.isdir(args.out_dir):
    os.makedirs(args.out_dir)

print("Options in effect:")
arguments = {}
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
    arguments[arg] = getattr(args, arg)
print("")

with open(os.path.join(args.out_dir, 'souporcell_settings.json'), 'w') as f:
    json.dump(arguments, f, indent=4)
f.close()


if args.no_umi == "True":
    args.no_umi = True
else:
    args.no_umi = False

print("checking modules")
# importing all reqs to make sure things are installed
import numpy as np
import scipy
import gzip
import math
import pystan
import vcf
import pysam
import pyfaidx
import subprocess
import time
import os
print("imports done")


open_function = lambda f: gzip.open(f,"rt") if f[-3:] == ".gz" else open(f)

print("checking bam for expected tags")
UMI_TAG = args.umi_tag
CELL_TAG = args.cell_tag
assert CELL_TAG == "CB", "vartrix doesnt support different cell tags, remake bam with cell tag as CB"
#load each file to make sure it is legit
bc_set = set()
with open_function(args.barcodes) as barcodes:
    for (index, line) in enumerate(barcodes):
        bc = line.strip()
        bc_set.add(bc)

assert len(bc_set) > 50, "Fewer than 50 barcodes in barcodes file? We expect 1 barcode per line."

assert not(not(args.known_genotypes == None) and not(args.common_variants == None)), "cannot set both know_genotypes and common_variants"
if args.known_genotypes_sample_names:
    assert not(args.known_genotypes == None), "if you specify known_genotype_sample_names, must specify known_genotypes option"
    assert len(args.known_genotypes_sample_names) == int(args.clusters), "length of known genotype sample names should be equal to k/clusters"
if args.known_genotypes:
    reader = vcf.Reader(open(args.known_genotypes))
    assert len(reader.samples) >= int(args.clusters), "number of samples in known genotype vcfs is less than k/clusters"
    if args.known_genotypes_sample_names == None:
        args.known_genotypes_sample_names = reader.samples
    for sample in args.known_genotypes_sample_names:
        assert sample in args.known_genotypes_sample_names, "not all samples in known genotype sample names option are in the known genotype samples vcf?"

#test bam load
bam = pysam.AlignmentFile(args.bam)
num_cb = 0
num_cb_cb = 0 # num reads with barcodes from barcodes.tsv file
num_umi = 0
num_read_test = 100000
for (index,read) in enumerate(bam):
    if index >= num_read_test:
        break
    if read.has_tag(CELL_TAG):
        num_cb += 1
        if read.get_tag(CELL_TAG) in bc_set:
            num_cb_cb += 1
    if read.has_tag(UMI_TAG):
        num_umi += 1
if args.skip_remap and args.common_variants == None and args.known_genotypes == None:
        assert False, "WARNING: skip_remap enables without common_variants or known genotypes. Variant calls will be of poorer quality. Turn on --ignore True to ignore this warning"

assert float(num_cb) / float(num_read_test) > 0.5, "Less than 50% of first 100000 reads have cell barcode tag (CB), turn on --ignore True to ignore"
if not(args.no_umi):
    assert float(num_umi) / float(num_read_test) > 0.5, "Less than 50% of first 100000 reads have UMI tag (UB), turn on --ignore True to ignore"
assert float(num_cb_cb) / float(num_read_test) > 0.05, "Less than 25% of first 100000 reads have cell barcodes from barcodes file, is this the correct barcode file? turn on --ignore True to ignore"

print("checking fasta")
fasta = pyfaidx.Fasta(args.fasta, key_function = lambda key: key.split()[0])

print("done")

#### END MAIN RUN SCRIPT
