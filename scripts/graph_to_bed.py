#!/usr/bin/env python3

import argparse
from collections import defaultdict
import os
import sys

from intervaltree import IntervalTree

if sys.version_info[0] < 3:
    raise Exception("Must be using Python3")


def read_graph(graphf):
    intD = defaultdict(IntervalTree)
    with open(graphf) as infile:
        for line in infile:
            if line.startswith("sequence"):
                fields = line.rstrip().rsplit()
                l, r = fields[1], fields[2]
                lchrom, lpos = l[:-1].rsplit(":")
                rchrom, rpos = r[:-1].rsplit(":")
                lpos, rpos = int(lpos), int(rpos) + 1  #use the +1 semi-closed coordinate system of UCSC Genome Browser.
                cn = float(fields[3])

                if add_chr_tag and not lchrom.startswith('chr'):
                    lchrom = "chr" + lchrom
                    rchrom = "chr" + rchrom

                intD[lchrom].addi(lpos,rpos,cn)

    return intD


def readFlist(filelist):
    flist = []
    with open(filelist) as infile:
        for line in infile:
            line = line.rstrip()
            if line:
                fields = line.rsplit()
                if len(fields) < 2 or len(fields) > 3:
                    print("Bad formatting in: ", line)
                else:
                    flist.append(fields)

    return flist

def merge_intervals(cn_segs):
    msegs = []

    if len(cn_segs) < 2:
        return cn_segs

    prev_interval = cn_segs[0]
    for i in range(1, len(cn_segs)):
        if prev_interval[0] != cn_segs[i][0] or cn_segs[i][1] - prev_interval[2] > 1:
            msegs.append(prev_interval)
            prev_interval = cn_segs[i]

        else:
            prev_interval[2] = cn_segs[i][2]

    msegs.append(prev_interval)

    return msegs


def make_bed(intD, min_cn = 0, unmerged = False):
    # sort by chromosome name
    sortablenames = []
    realnames = intD.keys()
    for x in realnames:
        stripName = x
        if x.startswith("chr"):
            stripName = x.lstrip("chr")

        try:
            stripName = int(stripName)

        except ValueError:
            pass

        sortablenames.append(stripName)

    orderedNames, _ = zip(*sorted(zip(realnames, sortablenames), key=lambda v: (isinstance(v[1], str), v[1])))

    cn_segs = []
    for c in orderedNames:
        curr_it = intD[c]
        ccs = sorted([[c, x.begin, x.end] for x in curr_it if x.data > min_cn])
        cn_segs.extend(ccs)

    bedlist = merge_intervals(cn_segs) if not unmerged else cn_segs
    return bedlist


def write_bed(bedlist, ofname):
    with open(ofname, 'w') as outfile:
        for t in bedlist:
            l = [str(x) for x in t]
            outfile.write("\t".join(l) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert segments in graph files to .bed format")
    parser.add_argument("-i", "--input", help="Path to list of files to use. Each line formatted as: \
    samplename /path/to/sample_amplicon1_cycles.txt /path/to/sample_amplicon1_graph.txt", required=True)
    parser.add_argument("--min_cn", type=float, help="Minimum CN to report region (default 0)", default=0)
    parser.add_argument("--unmerged", help="Do not merge adjacent intervals from graph file", action='store_true', default=False)
    parser.add_argument("--add_chr_tag", help="Add \'chr\' to the beginning of chromosome names in input files",
                        action='store_true', default=False)

    args = parser.parse_args()
    add_chr_tag = args.add_chr_tag

    collection_name = os.path.splitext(args.input)[0] + "_bed_files/"
    os.makedirs(collection_name, exist_ok=True)

    flist = readFlist(args.input)

    for fpair in flist:
        if len(fpair) > 2:
            sName, cyclesFile, graphFile = fpair

        else:
            print(fpair)
            sys.stderr.write("File list not properly formatted\n")
            sys.exit(1)

        ofname = collection_name + os.path.splitext(os.path.basename(graphFile))[0] + ".bed"
        print(ofname)

        intD = read_graph(graphFile)
        bedlist = make_bed(intD, args.min_cn, args.unmerged)
        write_bed(bedlist, ofname)
