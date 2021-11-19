#!/usr/bin/env python

import argparse
from collections import defaultdict
import copy
import os

from intervaltree import IntervalTree

# Break AA seeds on known low mappability, centromere, segmental duplication regions.
# This enables AA to better handle ultra-long AA seeds that can sometimes arise depending on the choice of CNV caller.


# read a bed file into a dictionary of interval trees, where keys are chromosomes
def read_bed(ifname, keepdat=False):
    beddict = defaultdict(IntervalTree)
    with open(ifname) as infile:
        for line in infile:
            line = line.rstrip()
            if line:
                fields = line.rsplit()
                s, e = int(fields[1]), int(fields[2])
                if e - s == 0:
                    print("Size 0 interval found. Skipping: " + line)
                    continue

                if keepdat:
                    beddict[fields[0]].addi(s, e, tuple(fields[3:]))
                else:
                    beddict[fields[0]].addi(s, e)

    return beddict


# write a bed file from a bed dict of interval trees, where keys are chromosomes. filters elements smaller than minsize
def write_bed(ofname, beddict, minsize):
    with open(ofname, 'w') as outfile:
        for k in sorted(beddict.keys()):
            ivallist = sorted([[k, x.begin, x.end] + list(x.data) for x in beddict[k]])
            for i in ivallist:
                if i[2] - i[1] >= minsize:
                    oline = "\t".join([str(x) for x in i]) + "\n"
                    outfile.write(oline)


# read regions to split on/filter into dictionary of interval trees, where keys are chromosomes
def read_filt_regions(ref):
    filt_regions = defaultdict(IntervalTree)
    AA_DATA_REPO = os.environ["AA_DATA_REPO"] + "/" + ref + "/"
    fdict = {}
    with open(AA_DATA_REPO + "file_list.txt") as infile:
        for line in infile:
            line = line.rstrip()
            if line:
                fields = line.rsplit()
                fdict[fields[0]] = fields[1]

    # do merged_centromeres_conserved
    mcc = AA_DATA_REPO + ref + "_merged_centromeres_conserved_sorted.bed"
    # do map excluded
    mex = AA_DATA_REPO + fdict["mapability_exclude_filename"]
    # do segdups
    sdp = AA_DATA_REPO + fdict["segdup_filename"]
    for cf in [mcc, mex, sdp]:
        print("Reading " + cf)
        result = read_bed(cf)
        for k, ivalt in result.items():
            filt_regions[k].update(ivalt)

    return filt_regions


# create a new seed tree of regions that have filtered regions trimmed out
def trim_seeds(seeddict, filt_regions):
    updated_seeds = defaultdict(IntervalTree)
    for k, ivalt in seeddict.items():
        updated_seed_tree = copy.copy(ivalt)
        for ival in ivalt:
            if ival.end - ival.begin > 1000000:
                for h in filt_regions[k][ival.begin: ival.end]:
                    updated_seed_tree.chop(h.begin - 250, h.end + 250)

        updated_seeds[k] = updated_seed_tree

    return updated_seeds


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Break AA seeds on elements which AA cannot analyze.")
    parser.add_argument("--ref", help="Reference genome version.", choices=["hg19", "GRCh37", "GRCh38", "mm10",
                        "GRCm38"], required=True)
    parser.add_argument("--seeds", help="path to bed file of seed regions", type=str, required=True)
    parser.add_argument("--minsize", help="Minimum trimmed seed size to keep", type=float, default=50000)
    args = parser.parse_args()

    p, f = os.path.split(args.seeds)
    outname = p + os.path.splitext(f)[0] + "_trimmed.bed"

    print("Reading seeds and filter regions")
    seeddict = read_bed(args.seeds, keepdat=True)
    filt_regions = read_filt_regions(args.ref)
    print("Updating seeds")
    updated_seeds = trim_seeds(seeddict, filt_regions)
    write_bed(outname, updated_seeds, args.minsize)
    print(outname)
