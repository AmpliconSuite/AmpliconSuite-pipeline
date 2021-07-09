#!/usr/bin/env python3

import argparse
from collections import defaultdict
import os
import sys

from intervaltree import IntervalTree


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


def buildregiond(regionlist):
    rd = defaultdict(IntervalTree)
    for r in regionlist:
        chrom, posp = r.rsplit(":")
        l, r = posp.rsplit("-")
        l, r = sorted((int(l), int(r)))
        rd[chrom].addi(l, r)

    return rd


def read_graph(graphf, rd, intD, sname):
    with open(graphf) as infile:
        for line in infile:
            if line.startswith("discordant"):
                fields = line.rstrip().rsplit()
                l,r = fields[1].rsplit("->")
                lchrom, lpos = l[:-1].rsplit(":")
                rchrom, rpos = r[:-1].rsplit(":")
                lpos, rpos = int(lpos), int(rpos)
                cn = float(fields[2])

                if add_chr_tag and not lchrom.startswith('chr'):
                    lchrom = "chr" + lchrom
                    rchrom = "chr" + rchrom

                hitA = len(rd[lchrom][lpos:lpos+1]) > 0
                hitB = len(rd[rchrom][rpos:rpos+1]) > 0

                extern = "external"
                if hitA or hitB:
                    if hitA and hitB:
                        extern = "internal"

                    spairs = sorted([(lchrom,lpos),(rchrom,rpos)])
                    intD[spairs[0][0]].append((spairs[0][0], spairs[0][1], spairs[1][0], spairs[1][1], extern, sname))
                    # intD[lchrom].append((lpos, rpos, extern))


def write_output(intD, outfile):
    if intD:
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
        for c in orderedNames:
            sorted_list = sorted(intD[c])
            for t in sorted_list:
                l = [str(x) for x in t]
                outfile.write("\t".join(l) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract breakpoints in graph file as a bed file")
    parser.add_argument("-i", "--input", help="Path to list of files to use. Each line formatted as: \
    samplename /path/to/sample_amplicon1_cycles.txt /path/to/sample_amplicon1_graph.txt", required=True)
    # parser.add_argument("--min_cn", type=float, help="Minimum CN to report region (default 0)", default=0)
    parser.add_argument("-r","--regions", help="restrict edge extraction to the following regions, formatted as \n"
                                               "-r chrA:start-stop chrB:start-stop ...", nargs="+", default=[])
    parser.add_argument("--add_chr_tag", help="Add \'chr\' to the beginning of chromosome names in input files",
                        action='store_true', default=False)

    args = parser.parse_args()
    add_chr_tag = args.add_chr_tag
    flist = readFlist(args.input)
    regionstring = "_".join(args.regions) if args.regions else "all"
    ofname = os.path.splitext(os.path.basename(args.input))[0] + "_" + regionstring + "_breakpoints.bed"
    ofname = ofname.replace(':', '-')
    print(ofname)
    with open(ofname, 'w') as outfile:
        intD = defaultdict(list)

        for fpair in flist:
            if len(fpair) > 2:
                sName, cyclesFile, graphFile = fpair

            else:
                print(fpair)
                sys.stderr.write("File list not properly formatted\n")
                sys.exit(1)

            rd = buildregiond(args.regions)
            sname = os.path.splitext(os.path.basename(graphFile))[0].rsplit("_graph")[0]
            read_graph(graphFile, rd, intD, sname)

        write_output(intD, outfile)
