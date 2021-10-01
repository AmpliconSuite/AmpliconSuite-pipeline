#!/usr/bin/env python3

# author: Jens Luebeck (jluebeck [at] ucsd.edu)

import argparse
from collections import defaultdict
import os
import sys

from intervaltree import IntervalTree

if sys.version_info[0] < 3:
    raise Exception("Must be using Python3")


def read_graph(graphf):
    intD = defaultdict(IntervalTree)
    de_list = []
    with open(graphf) as infile:
        for line in infile:
            if line.startswith("sequence"):
                fields = line.rstrip().rsplit()
                l, r = fields[1], fields[2]
                lchrom, lpos = l[:-1].rsplit(":")
                rchrom, rpos = r[:-1].rsplit(":")
                lpos, rpos = int(lpos), int(rpos)
                if lpos == rpos: rpos += 1
                cn = float(fields[3])

                if add_chr_tag and not lchrom.startswith('chr'):
                    lchrom = "chr" + lchrom
                    rchrom = "chr" + rchrom

                intD[lchrom].addi(lpos, rpos, cn)

            elif line.startswith("discordant"):
                fields = line.rstrip().rsplit()
                l, r = fields[1].rsplit("->")
                lchrom, lpos = l[:-1].rsplit(":")
                rchrom, rpos = r[:-1].rsplit(":")
                lpos, rpos = int(lpos), int(rpos)
                strand1 = l[-1]
                strand2 = '-' if r[-1] == '+' else '+'
                if add_chr_tag and not lchrom.startswith('chr'):
                    lchrom = "chr" + lchrom
                    rchrom = "chr" + rchrom

                de_list.append((lchrom, lpos, rchrom, rpos, strand1, strand2))

    return intD, de_list


def readFlist(filelist):
    flist = []
    with open(filelist) as infile:
        for line in infile:
            line = line.rstrip()
            if line:
                fields = line.rsplit()
                flist.append((fields[0], fields[-1]))

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


def make_bed(intD, min_cn=0, unmerged=False):
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
        if unmerged:
            ccs = sorted([[c, x.begin, x.end, x.data] for x in curr_it if x.data > min_cn])
        else:
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
    parser = argparse.ArgumentParser(description="Convert segments in graph files to .bed format. If CN column is "
                                                 "desired, must set --unmerged.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-g", "--graph", help="Path to graph file. Either --graph or --input specifying a list are "
                                             "required", type=str)
    group.add_argument("-i", "--input", help="Path to list of files to use. Each line containing the path to a graph "
                                             "file in the last column", type=str)
    parser.add_argument("--min_cn", type=float, help="Minimum CN to report region (default 0)", default=0.0)
    parser.add_argument("--unmerged", help="Do not merge adjacent intervals from graph file", action='store_true',
                        default=False)
    parser.add_argument("--add_chr_tag", help="Add \'chr\' to the beginning of chromosome names in input files",
                        action='store_true', default=False)
    # parser.add_argument("--graph_edge_bedpe", help="Also report a .bedpe file of all breakpoint graph edges (default "
    #                                                "true)", action='store_true', default=True)

    args = parser.parse_args()
    add_chr_tag = args.add_chr_tag

    if args.input:
        collection_name = os.path.splitext(args.input)[0] + "_bed_files/"
        os.makedirs(collection_name, exist_ok=True)
        flist = readFlist(args.input)

    else:
        collection_name = ""
        sname = os.path.splitext(os.path.basename(args.graph))[0]
        flist = [(sname, args.graph)]

    for fpair in flist:
        sName, graphFile = fpair
        intD, de_list = read_graph(graphFile)
        bedlist = make_bed(intD, args.min_cn, args.unmerged)
        ofname = collection_name + os.path.splitext(os.path.basename(graphFile))[0] + ".bed"
        print(ofname)
        write_bed(bedlist, ofname)
        ofname = collection_name + os.path.splitext(os.path.basename(graphFile))[0] + "_breakpoints.bedpe"
        print(ofname)
        write_bed(de_list, ofname)
