#!/usr/bin/env python

import argparse
from collections import defaultdict
import os
import sys

from intervaltree import Interval, IntervalTree

h1 = "SequenceEdge: StartPosition, EndPosition, PredictedCopyCount, AverageCoverage, Size, NumberReadsMapped\n"
h2 = "BreakpointEdge: StartPosition->EndPosition, PredictedCopyCount, NumberOfReadPairs, " \
     "HomologySizeIfAvailable(<0ForInsertions), Homology/InsertionSequence\n"
max_hop = 800


# read graph file
def read_graph(graphf, maxhopsize, filter_non_everted, max_support):
    intD = defaultdict(IntervalTree)
    edge_line_list = []
    dbp_set = set()
    removed_count = 0
    with open(graphf) as infile:
        for line in infile:
            if line.startswith("sequence"):
                fields = line.rstrip().rsplit()
                l, r = fields[1], fields[2]
                lchrom, lpos = l[:-1].rsplit(":")
                rchrom, rpos = r[:-1].rsplit(":")
                if int(lpos) - int(rpos) != 0:
                    lpos, rpos = int(lpos), int(rpos)
                else:
                    lpos, rpos = int(lpos), int(rpos) + 1

                cn = float(fields[3])
                cov = float(fields[4])
                seg_size = int(fields[5])
                nrm = int(fields[6])
                intD[lchrom].addi(lpos, rpos, (cn, cov, seg_size, nrm))

            elif any([line.startswith(x) for x in ["source", "concordant", "discordant"]]):
                fields = line.rstrip().rsplit("\t")
                pair = fields[1].rsplit("->")
                left = pair[0].rsplit(":")
                lchrom, lpos, ldir = left[0], int(left[1][:-1]), left[1][-1]
                right = pair[1].rsplit(":")
                rchrom, rpos, rdir = right[0], int(right[1][:-1]), right[1][-1]
                support = 1
                if fields[0] == "discordant":
                    support = int(fields[3])

                # expected orientation: ldir == "+" and rdir == "-", on concordant but backwards for discordant:
                if support <= max_support:
                    if fields[0] == "discordant" and rchrom == lchrom and 0 < lpos - rpos <= maxhopsize and \
                            ldir == '+' and rdir == '-':
                        print("Removing: " + line.rstrip() + " | hopsize: " + str(abs(lpos - rpos)))
                        removed_count+=1
                        continue

                    if fields[0] == "discordant" and rchrom == lchrom and 0 < lpos - rpos <= maxhopsize and \
                            ldir == '-' and rdir == '+' and filter_non_everted:
                        print("Removing: " + line.rstrip() + " | hopsize: " + str(abs(lpos - rpos)))
                        removed_count+=1
                        continue

                edge_line_list.append(line)
                if line.startswith("discordant"):
                    dbp_set.add((lchrom, lpos))
                    dbp_set.add((rchrom, rpos))

    return intD, edge_line_list, dbp_set, removed_count


def write_graph(outname, merged_seg_intd, edge_line_list):
    sorder_chroms = ["chr" + str(x) for x in range(1, 23)] + ["chrX", "chrY"] + [str(x) for x in range(1, 23)] + ["X,Y"]
    with open(outname, 'w') as outfile:
        outfile.write(h1)
        for chrom in sorder_chroms:
            svals = merged_seg_intd[chrom]
            sorted_svals = sorted(svals, key=lambda x: (x.begin, x.end, x.data))
            for x in sorted_svals:
                outline = "\t".join([str(y) for y in ["sequence", chrom + ":" + str(x[0]) + "-", chrom + ":" + str(x[1]) + "+"] + list(x[2])]) + "\n"
                outfile.write(outline)

        outfile.write(h2)
        for eline in edge_line_list:
            outfile.write(eline)


# return true if both ends clear of breakpoints and it's short
def is_orphan(ival, chrom, dbp_set):
    has_bps = ((chrom, ival.begin) in dbp_set or (chrom, ival.end) in dbp_set)
    return not has_bps and ival.end - ival.begin <= max_hop


# compute proportion of graph segments over size cutoff
def proportion_over_size(seg_intd, size=25000):
    total = 0.0
    n_over = 0.0
    for chrom, intervals in seg_intd.items():
        for x in intervals:
            s = x.end - x.begin
            if s > size:
                n_over += 1
            total += 1

    return n_over/(max(1.0, total))


def ClusterIntervalsFromSortedList(seg_intd, dbp_set):
    # Sorting based on the increasing order
    # of the start intervals
    clustered_ivald = defaultdict(list)
    for chrom, unsort_intervals in seg_intd.items():
        intervals = sorted(unsort_intervals, key=lambda x: x.begin)
        if len(intervals) < 2:
            clustered_ivald[chrom] = [intervals]

        clust_list = []
        #set an initial
        a = intervals[0]
        curr_clust = [a]
        for b in intervals[1:]:
            can_merge = False
            has_orphan = is_orphan(a, chrom, dbp_set) or is_orphan(b, chrom, dbp_set)
            if b.begin - a.end == 1 and has_orphan:
                if not (chrom, a.end) in dbp_set and not (chrom, b.begin) in dbp_set:  # endpoints can be joined
                    can_merge = True
                    curr_clust.append(b)

            if not can_merge:
                clust_list.append(curr_clust)
                curr_clust = [b]

            a = b

        clust_list.append(curr_clust)
        # debug printing
        # for x in clust_list:
        #     if len(x) < 2:
        #         continue
        #     for y in x:
        #         sys.stdout.write(str(y.begin) + "-" + str(y.end) + "|" + str(y.end - y.begin) + "|" +
        #                          str(is_orphan(y, chrom, dbp_set)) + "\t")
        #
        #     sys.stdout.write("\n")

        clustered_ivald[chrom] = clust_list

    return clustered_ivald


def merge_clusters(clustered_ivald):
    merged_seg_intd = defaultdict(list)
    for chrom, cluster_list in clustered_ivald.items():
        print(chrom, "merging orphans")
        merged_seg_list = []
        for c in cluster_list:
            #identify the orphaned segments, identify the parents. merge appropriately
            #need to handle case where two different CNs present in parents
            is_parent_vect = [False]*len(c)
            for ind, x in enumerate(c):
                if x.end - x.begin > max_hop:
                    is_parent_vect[ind] = True

            curr_orphan_chain = []
            orphan_collapsed_intervals = []
            collapsed_is_parent = []
            for ip, x in zip(is_parent_vect, c):
                if ip:
                    if curr_orphan_chain:
                        # merge
                        lens = [i.end - i.begin for i in curr_orphan_chain]
                        cns = [j.data[0] for j in curr_orphan_chain]
                        covs = [j.data[1] for j in curr_orphan_chain]
                        tsize = sum([j.data[2] for j in curr_orphan_chain])
                        tnrm = sum([j.data[3] for j in curr_orphan_chain])
                        total_cn_weight = sum([i*j for i, j in zip(lens, cns)])
                        total_cov_weight = sum([i*j for i, j in zip(lens, covs)])
                        mean_cn_weight = total_cn_weight/sum(lens)
                        mean_cov_weight = total_cov_weight/sum(lens)
                        nd = (mean_cn_weight, mean_cov_weight, tsize, tnrm)
                        nn = Interval(curr_orphan_chain[0].begin, curr_orphan_chain[-1].end, nd)
                        orphan_collapsed_intervals.append(nn)
                        collapsed_is_parent.append(False)
                        curr_orphan_chain = []

                    orphan_collapsed_intervals.append(x)
                    collapsed_is_parent.append(True)

                else:
                    # put orphan on chain
                    curr_orphan_chain.append(x)

            # fencepost
            if curr_orphan_chain:
                # merge
                lens = [i.end - i.begin for i in curr_orphan_chain]
                cns = [j.data[0] for j in curr_orphan_chain]
                covs = [j.data[1] for j in curr_orphan_chain]
                tsize = sum([j.data[2] for j in curr_orphan_chain])
                tnrm = sum([j.data[3] for j in curr_orphan_chain])
                total_cn_weight = sum([i * j for i, j in zip(lens, cns)])
                total_cov_weight = sum([i * j for i, j in zip(lens, covs)])
                mean_cn_weight = total_cn_weight / sum(lens)
                mean_cov_weight = total_cov_weight / sum(lens)
                nd = (mean_cn_weight, mean_cov_weight, tsize, tnrm)
                nn = Interval(curr_orphan_chain[0].begin, curr_orphan_chain[-1].end, nd)
                orphan_collapsed_intervals.append(nn)
                collapsed_is_parent.append(False)

            # debug printing
            if len(orphan_collapsed_intervals) > 1:
                for yind, y in enumerate(orphan_collapsed_intervals):
                    sys.stdout.write(str(y.begin) + "-" + str(y.end) + "|" + str(y.end - y.begin) + "|" +
                                     str(not collapsed_is_parent[yind]) + "|" + str(y.data[0]) + "\t")

                sys.stdout.write("\n")

            # Now assign orphan to parents. Join parents if needed.
            join_to_next = [False]*len(orphan_collapsed_intervals)
            if len(orphan_collapsed_intervals) < 2:
                fully_joined = orphan_collapsed_intervals

            else:
                ind = 0
                # print("final merging")
                for cip, x in zip(collapsed_is_parent, orphan_collapsed_intervals):
                    # print(cip,x)
                    join_both_parents = False
                    join_dir = None
                    if not cip:
                        if 0 < ind < len(collapsed_is_parent) - 1:  # it's a middle
                            # print("MIDDLE")
                            if abs(orphan_collapsed_intervals[ind-1].data[0] - orphan_collapsed_intervals[ind+1].data[0]) < 1:
                                join_both_parents = True
                                join_dir = 0

                            else:
                                if abs(orphan_collapsed_intervals[ind-1].data[0] - x.data[0]) <= abs(orphan_collapsed_intervals[ind+1].data[0] - x.data[0]):
                                    join_dir = -1
                                else:
                                    join_dir = 0

                        elif ind == 0:
                            join_dir = 0

                        elif ind == len(collapsed_is_parent) - 1:
                            join_dir = -1

                    # print("JD", join_dir)
                    if join_both_parents:
                        hits = [-1, 0]

                    elif join_dir is not None:
                        hits = [join_dir]

                    else:
                        # print("EMPTY HITS")
                        hits = []

                    for h in hits:
                        join_to_next[ind + h] = True

                    ind += 1
                    # print(ind,"J2N",join_to_next)

                # print(join_to_next)
                fully_joined = []
                curr_joins = []
                for join, ival in zip(join_to_next, orphan_collapsed_intervals):
                    if not join:
                        # if curr_joins:
                        # merge
                        curr_joins.append(ival)
                        lens = [i.end - i.begin for i in curr_joins]
                        cns = [j.data[0] for j in curr_joins]
                        covs = [j.data[1] for j in curr_joins]
                        tsize = sum([j.data[2] for j in curr_joins])
                        tnrm = sum([j.data[3] for j in curr_joins])
                        total_cn_weight = sum([i * j for i, j in zip(lens, cns)])
                        total_cov_weight = sum([i * j for i, j in zip(lens, covs)])
                        mean_cn_weight = total_cn_weight / sum(lens)
                        mean_cov_weight = total_cov_weight / sum(lens)
                        nd = (mean_cn_weight, mean_cov_weight, tsize, tnrm)
                        nn = Interval(curr_joins[0].begin, curr_joins[-1].end, nd)
                        fully_joined.append(nn)
                        curr_joins = []

                        # fully_joined.append(ival)

                    else:
                        # put orphan on chain
                        curr_joins.append(ival)

                # fencepost
                if curr_joins:
                    # merge
                    lens = [i.end - i.begin for i in curr_joins]
                    cns = [j.data[0] for j in curr_joins]
                    covs = [j.data[1] for j in curr_joins]
                    tsize = sum([j.data[2] for j in curr_joins])
                    tnrm = sum([j.data[3] for j in curr_joins])
                    total_cn_weight = sum([i * j for i, j in zip(lens, cns)])
                    total_cov_weight = sum([i * j for i, j in zip(lens, covs)])
                    mean_cn_weight = total_cn_weight / sum(lens)
                    mean_cov_weight = total_cov_weight / sum(lens)
                    nd = (mean_cn_weight, mean_cov_weight, tsize, tnrm)
                    nn = Interval(curr_joins[0].begin, curr_joins[-1].end, nd)
                    fully_joined.append(nn)

            merged_seg_list.extend(fully_joined)
            # debug printing
            if len(orphan_collapsed_intervals) > 1:
                print("fully joined")
                for yind, y in enumerate(fully_joined):
                    sys.stdout.write(str(y.begin) + "-" + str(y.end) + "|" + str(y.end - y.begin) + "|" + str(y.data[0]) +  "\t")

                sys.stdout.write("\n")
                sys.stdout.write("\n")

        merged_seg_intd[chrom] = merged_seg_list

    return merged_seg_intd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Remove spurious sequence artifact edges from AA graph")
    parser.add_argument("-g", help="path to graph", type=str)
    parser.add_argument("--graph_list", help="Text file listing paths to all graphs to clean", type=str)
    parser.add_argument("--max_hop_size", help="Maximum size of everted read hop (default 4000)", type=float,
                        default=4000)
    parser.add_argument("--max_hop_support", help="Maximum number of discordant read pairs to consider as hop "
                                                  "(default (10)", type=int, default=10)
    parser.add_argument("--filter_non_everted", help="Filter non-everted hops", action='store_true', default=False)

    args = parser.parse_args()

    if not args.g and not args.graph_list:
        sys.stderr.write("At least one of -g or --graph_list required.\n")
        sys.exit(1)

    elif args.g:
        gfile_list = [args.g]

    else:
        with open(args.graph_list) as f:
            gfile_list = f.read().splitlines()

    max_hop = args.max_hop_size

    for gfile in gfile_list:
        print("Cleaning " + gfile)
        p, f = os.path.split(gfile)
        outname = f.rsplit("_graph.txt")[0] + "_cleaned_graph.txt"
        seg_intd, edge_line_list, dbp_set, removed_count = read_graph(gfile, max_hop, args.filter_non_everted,
                                                                      args.max_hop_support)
        # compute the fraction of segments over 50 kbp
        print("Initial proportion over 25kbp:", proportion_over_size(seg_intd))
        if removed_count > 0:
            clustered_ivald = ClusterIntervalsFromSortedList(seg_intd, dbp_set)
            merged_ivald = merge_clusters(clustered_ivald)
            print("Merged proportion over 25kbp", proportion_over_size(merged_ivald))

        else:
            print("Nothing to clean")
            merged_ivald = seg_intd

        write_graph(outname, merged_ivald, edge_line_list)
        print(outname + "\n")
