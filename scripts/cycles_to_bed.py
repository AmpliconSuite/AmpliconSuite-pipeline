#!/usr/bin/env python

import argparse
from collections import defaultdict
import os

from intervaltree import IntervalTree

# Write cycles files into merged beds

# read a cycles file into a dictionary of interval tree dictionaries
def read_cycles_file(fname):
    # all_cycles maps cycle_num -> chromosome -> IntervalTree
    seglookup = {}
    all_cycles_ivald = {}
    all_cycles = {}
    with open(fname) as infile:
        for line in infile:
            if line.startswith("Segment"):
                fields = line.rstrip().rsplit()
                segnum = int(fields[1])
                dtup = (fields[2], int(fields[3]), int(fields[4]))
                seglookup[segnum] = dtup

            elif line.startswith(("Cycle")):
                cycle_ivalt = defaultdict(IntervalTree)
                slist = []
                fields = line.rstrip().rsplit(';')
                cd = {x.rsplit("=")[0]: x.rsplit("=")[1] for x in fields}
                cnum = int(cd["Cycle"])
                segstring = cd["Segments"]
                copy_count = cd["Copy_count"]
                seglist = [(int(x[:-1]), x[-1]) for x in segstring.rsplit(",")]
                for x in seglist:
                    if x[0] == 0:
                        continue

                    dtup = seglookup[x[0]]
                    slist.append(dtup + (x[1], copy_count))
                    cycle_ivalt[dtup[0]].addi(dtup[1], dtup[2], (x[1], copy_count))  # duplicates will be overwritten

                all_cycles[cnum] = slist
                all_cycles_ivald[cnum] = cycle_ivalt

    return all_cycles_ivald, all_cycles


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
            if a[0] > max + 1:
                if i != 0:
                    m.append([chrom_num, s, max])
                max = a[1]
                s = a[0]
            else:
                if a[1] > max:
                    max = a[1]

        # 'max' value gives the last point of
        # that particular interval
        # 's' gives the starting point of that interval
        # 'm' array contains the list of all merged intervals
        if max != -100000 and [chrom_num, s, max] not in m:
            m.append([chrom_num, s, max])

        merge_ivals.extend(m)

    s_ivals = sorted(merge_ivals, key=lambda x: x[0])
    for i in range(len(s_ivals)):
        s_ivals[i][0] = 'chr' + str(s_ivals[i][0])

    return s_ivals


def write_bed(prefix, merged_intervals):
    with open(prefix + ".bed", 'w') as outfile:
        for i in merged_intervals:
            oline = "\t".join([str(x) for x in i]) + "\n"
            outfile.write(oline)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert cycles file to a series of bed files")
    parser.add_argument("-c", "--cycles", help="Cycles file", type=str, required=True)
    args = parser.parse_args()

    a_cycles_ivald, a_cycles_list = read_cycles_file(args.cycles)
    fpre = os.path.splitext(os.path.basename(args.cycles))[0]
    for cnum, curr_cycle in a_cycles_ivald.items():
        pre = fpre + "_unordered_cycle" + str(cnum)
        merged_intervals = mergeIntervals(curr_cycle)
        write_bed(pre, merged_intervals)
        pre = fpre + "_cycle" + str(cnum)
        write_bed(pre, a_cycles_list[cnum])
