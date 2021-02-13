#!/usr/bin/env python

import argparse
from collections import defaultdict
import os

from intervaltree import IntervalTree

# Write cycles files into merged beds

# read a cycles file into a a dictionary of interval tree dictionaries
def read_cycles_file(fname):
    # all_cycles maps cycle_num -> chromosome -> IntervalTree
    seglookup = {}
    all_cycles = {}
    with open(fname) as infile:
        for line in infile:
            if line.startswith("Segment"):
                fields = line.rstrip().rsplit()
                segnum = int(fields[1])
                dtup = (fields[2], float(fields[3]), float(fields[4]))
                seglookup[segnum] = dtup

            elif line.startswith(("Cycle")):
                cycle_ivalt = defaultdict(IntervalTree)
                fields = line.rstrip().rsplit(';')
                cd = {x.rsplit("=")[0]: x.rsplit("=")[1] for x in fields}
                cnum = int(cd["Cycle"])
                segstring = cd["Segments"]
                seglist = [int(x[:-1]) for x in segstring.rsplit(",")]
                for x in seglist:
                    if x == 0:
                        continue

                    dtup = seglookup[x]
                    cycle_ivalt[dtup[0]].addi(dtup[1], dtup[2])  # duplicates will be overwritten

                all_cycles[cnum] = cycle_ivalt

    return all_cycles


def mergeIntervals(curr_cycle):
    merge_ivals = []

    for chrom, ivalt in curr_cycle.items():
        try:
            chrom_num = int(chrom.lstrip('chr'))
        except ValueError:
            chrom_num = chrom.lstrip('chr')

        # Sorting based on the increasing order
        # of the start intervals
        arr = sorted([(x.begin, x.end) for x in ivalt])
        # array to hold the merged intervals
        m = []
        s = -10000
        max = -100000
        for i in range(len(arr)):
            a = arr[i]
            if a[0] > max:
                if i != 0:
                    m.append([chrom_num, s, max])
                max = a[1]
                s = a[0]
            else:
                if a[1] >= max:
                    max = a[1]

        # 'max' value gives the last point of
        # that particular interval
        # 's' gives the starting point of that interval
        # 'm' array contains the list of all merged intervals
        if max != -100000 and [chrom_num, s, max] not in m:
            m.append([chrom_num, s, max])

        merge_ivals.extend(m)

    return sorted(merge_ivals, key=lambda x: x[0])


def write_bed(prefix, merged_intervals):
    with open(prefix + ".bed", 'w') as outfile:
        for i in merged_intervals:
            oline = "\t".join(['chr' + str(i[0]), str(int(i[1])), str(int(i[2]))]) + "\n"
            outfile.write(oline)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert cycles file to a series of bed files")
    parser.add_argument("-c", "--cycles", help="Cycles file", type=str, required=True)
    args = parser.parse_args()

    a_cycles = read_cycles_file(args.cycles)
    fpre = os.path.splitext(os.path.basename(args.cycles))[0]
    for cnum, curr_cycle in a_cycles.items():
        pre = fpre + "_cycle" + str(cnum)
        merged_intervals = mergeIntervals(curr_cycle)
        write_bed(pre, merged_intervals)
