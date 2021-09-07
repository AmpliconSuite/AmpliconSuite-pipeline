#!/usr/bin/env python3

# author: Jens Luebeck (jluebeck [at] ucsd.edu)

import os
import copy
import argparse
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np


id_to_coords = {}
end_to_id = {}
start_to_id = {}

# edge set for each node (oriented)
edgeDict = defaultdict(set)

# segment raw copy numbers
raw_cn = defaultdict(float)
len_zero_segs = set()


# DFS recursion
def DFSUtil(v, currPath, usedCN, lcp):
    # Mark the current node as visited
    usedCN[abs(v)] += 1
    currPath.append(v)

    # Recur for all the vertices
    # adjacent to this vertex
    cLPath = currPath
    currEdgeSet = edgeDict[v]
    if currPath[0] in edgeDict[v]:
        lcp = currPath

    for i in currEdgeSet:
        if usedCN[abs(i)] < scaled_cns[abs(i)]:
            retPath, clcp = DFSUtil(i, copy.copy(currPath), copy.copy(usedCN), copy.copy(lcp))
            if len(retPath) > len(cLPath):
                cLPath = retPath

            if len(clcp) > len(lcp):
                lcp = clcp

    return cLPath, lcp


# The function to do DFS traversal. It uses recursive DFSUtil() 
def DFS(v):
    usedCN = [0] * len(edgeDict)
    currPath, lcp = [], []

    # Call the recursive helper function
    lp, lcp = DFSUtil(v, currPath, usedCN, lcp)
    return lp, lcp


def remove_duplicate_paths(candidates):
    kept = []
    for xf in candidates:
        xr = [-x for x in xf[::-1]]
        xfstr = "".join([str(x) for x in xf])
        xrstr = "".join([str(x) for x in xr])
        keepable = True
        for yh in kept:
            yfull = "".join([str(y) for y in yh + yh])
            if xfstr in yfull or xrstr in yfull:
                keepable = False
                break

        if keepable:
            kept.append(xf)

    return kept


# read the graph and make dictionary of edges
def read_graph(graphf, skip_short_jumps):
    with open(graphf) as infile:
        seqN = 0
        for line in infile:
            if line.startswith("sequence"):
                seqN += 1
                fields = line.rstrip().rsplit("\t")
                chrom = fields[1].rsplit(":")[0]
                p1 = int(fields[1].rsplit(":")[1][:-1])
                p2 = int(fields[2].rsplit(":")[1][:-1])
                raw_cn[seqN] = float(fields[3])

                id_to_coords[seqN] = (chrom, p1, p2)
                id_to_coords[-1 * seqN] = (chrom, p2, p1)

                end_to_id[(chrom, p2)] = seqN
                start_to_id[(chrom, p1)] = seqN

                if p1 == p2:
                    print("zero len segment", line)
                    len_zero_segs.add(seqN)

                else:
                    end_to_id[(chrom, p1)] = -1 * seqN
                    start_to_id[(chrom, p2)] = -1 * seqN

            elif line.startswith("concordant") or line.startswith("discordant"):
                fields = line.rstrip().rsplit("\t")
                pair = fields[1].rsplit("->")
                left = pair[0].rsplit(":")
                lchrom, lpos, ldir = left[0], int(left[1][:-1]), left[1][-1]

                right = pair[1].rsplit(":")
                rchrom, rpos, rdir = right[0], int(right[1][:-1]), right[1][-1]
                try:
                    leftI, rightI = end_to_id[(lchrom, lpos)], start_to_id[(rchrom, rpos)]
                except KeyError:
                    print("WARNING", lchrom, lpos, rchrom, rpos, "edge doesn't link known segments")
                    continue

                # expected orientation
                # if ldir == "+" and rdir == "-":
                if skip_short_jumps and fields[0] == "discordant" and rchrom == lchrom and abs(
                        rpos - lpos) < 800 and rdir != ldir:
                    print("removing", line.rstrip())
                    continue

                edgeDict[leftI].add(rightI)
                edgeDict[-1 * rightI].add(-1 * leftI)

                if leftI in len_zero_segs:
                    edgeDict[-1 * leftI].add(rightI)
                    edgeDict[rightI].add(-1 * leftI)

                elif rightI in len_zero_segs:
                    edgeDict[leftI].add(-1 * rightI)
                    edgeDict[-1 * rightI].add(leftI)

            # print(fields[1], leftI, rightI)


def get_scaled_cns(raw_cn, scaling_factor):
    scaled_cns = {}
    for s, c in raw_cn.items():
        ratio = c / scaling_factor
        if ratio < 0.5:
            scaled_cns[s] = 0
            print("dropped", s, c)

        else:
            scaled_cns[s] = int(np.round(ratio))
        # print(s, int(np.round(ratio)), c)

    return scaled_cns


def get_median_cn(min_cn_cutoff, runmode, min_cn_seg_size=100):
    useable_cns = [y for x,y in raw_cn.items() if abs(id_to_coords[x][2] - id_to_coords[x][1]) > min_cn_seg_size]
    m_u_cn = max(useable_cns)
    if runmode == 'isolated':
        p10 = m_u_cn/10.0 if m_u_cn > 10000 else m_u_cn/20.0
        min_cn_cutoff = max(p10, min_cn_cutoff)

    print("Setting min cn cutoff to " + str(min_cn_cutoff))
    #return np.median([x for x in useable_cns if x > min_cn_cutoff])
    return np.percentile([x for x in useable_cns if x > min_cn_cutoff], 40)


def write_cycles_file(paths, id_to_coords, pweights, scaling_factor, ofname, plens):
    with open(ofname, 'w') as outfile:
        postups = [v for k, v in sorted(id_to_coords.items()) if k > 0]
        pchrom, pstart, pend, iind, sind = postups[0][0], postups[0][1], postups[0][2], 1, 1
        seglines = [["Segment", str(sind), pchrom, str(pstart), str(pend)]]
        for postup in postups[1:]:
            sind += 1
            seglines.append(["Segment", str(sind)] + [str(x) for x in postup])
            if postup[0] != pchrom or postup[1] - pend > 1:
                olist = ["Interval", str(iind), pchrom, str(pstart), str(pend)]
                outfile.write("\t".join(olist) + "\n")
                pchrom, pstart, pend = postup
                iind += 1

            else:
                pend = postup[2]

        olist = ["Interval", str(iind), pchrom, str(pstart), str(pend)]
        outfile.write("\t".join(olist) + "\n")
        outfile.write("List of cycle segments\n")
        for sl in seglines:
            outfile.write("\t".join(sl) + "\n")

        for pind, p in enumerate(paths):
            if p:
                fmtP = [str(abs(x)) + "+" if x > 0 else str(abs(x)) + "-" for x in p]
                if p[0] not in edgeDict[p[-1]]:
                    fmtP = ["0+", ] + fmtP + ["0-", ]

                outfile.write(
                    "Cycle=" + str(pind + 1) + ";Copy_count=" + str(scaling_factor) + ";ProportionAmplifiedExplained=" +
                    pweights[pind] + ";Segments=" + ",".join(fmtP) + ";Length=" + str(plens[pind]) + "bp" + "\n")


def get_cmap(n, name='jet'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)


def plot_cn_and_multiplicity():
    us_id, us_scn, us_rcn, us_cols = [], [], [], []
    cmap = get_cmap(len(set(scaled_cns.values())))
    for x,y in scaled_cns.items():
        us_id.append(int(x))
        us_scn.append(y)
        us_rcn.append(raw_cn[x])
        us_cols.append(cmap(y))

    sorted_cns, sorted_ids, sorted_mults, sorted_cols = zip(*sorted(zip(us_rcn, us_id, us_scn, us_cols)))
    # plt.style.use('ggplot')
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    x_pos = list(range(len(sorted_ids)))

    ax1.bar(x_pos, sorted_cns, color=sorted_cols)

    ax1.set_ylim(bottom=0)
    ymin, ymax = ax1.get_ylim()
    ax2_yticks = [x*scaling_factor for x in range(0, int(max(sorted_mults)) + 1)]
    ax1.set_ylim(bottom=0, top = max(ax2_yticks[-1], ymax))
    ymin, ymax = ax1.get_ylim()
    ax2.set_ylim(bottom=0, top=ymax)
    ax2.set_yticks(ax2_yticks)
    ax2.set_yticklabels([str(x) for x in range(0, int(max(sorted_mults)) + 1)])

    for x, m in zip(x_pos, sorted_mults):
        ax1.hlines(m*scaling_factor, x, x+1, color='k')
        ax1.hlines(m*scaling_factor - 0.5*scaling_factor, x, x+1, color='lightgrey',linestyles='dashed')
        ax1.hlines(m*scaling_factor + 0.5*scaling_factor, x, x+1, color='lightgrey',linestyles='dashed')

    plt.xlabel("Segment ID")
    plt.xticks(x_pos, sorted_ids)
    # plt.xticks(rotation=90, fontsize=1)
    ax1.set_xticklabels(sorted_ids, rotation=90, fontsize=3)
    plt.savefig(ofpre + "_CN_multiplicity_plot.png", dpi=300)


parser = argparse.ArgumentParser(
    description="Attempt to identify a longest path in the breakpoint graph consistent with CN ratios")
parser.add_argument("-g", "--graph", help="AA-formatted graph file", required=True)
parser.add_argument("--scaling_factor", help="Estimated CN of elements appearing once in the amplicon, defaults to "
                                             "median CN in graph however highly recommend user provide this value based"
                                             " on manual inspection of the graph", type=float)
parser.add_argument("--remove_short_jumps",
                    help="Remove very short discordant edges ( < 800 bp) which are not inversions", action="store_true")
parser.add_argument("--keep_all_LC", help="Keep all longest cyclic paths of same score", action="store_true",
                    default=False)
parser.add_argument("--minimum_cn_for_median_calculation", help="If not setting scaling factor manually, set a minimum "
                                                                "copy number for the scaling factor median amplified CN"
                                                                "calculation (default 4.5)", type=float, default=4.5)
# parser.add_argument("--length_estimate", help="Target size of the path (in kbp)", type=float)
parser.add_argument("--runmode", help="Default mode is 'bulk' for bulk WGS. Also supports 'isolated' mode for PFGE or "
                        "targeted sequencing.", choices=['bulk', 'isolated'])

args = parser.parse_args()

read_graph(args.graph, args.remove_short_jumps)
# print(len(edgeDict)/2, "edges will be considered")

ofpre = os.path.basename(args.graph).rsplit("_graph.txt")[0]
ofname = ofpre + "_candidate_cycles.txt"

if args.scaling_factor:
    scaling_factor = args.scaling_factor

else:
    print("using median as scaling factor")
    # scaling_factor = np.median([x for x in raw_cn.values() if x > args.minimum_cn_for_median_calculation])
    scaling_factor = get_median_cn(args.minimum_cn_for_median_calculation, args.runmode)

print("scaling factor: ", scaling_factor)
scaled_cns = get_scaled_cns(raw_cn, scaling_factor)

# An alternate way to do the exploration, using only the node with highest degree as starting position
# maxSeg = None
# maxEC = 0
# for s, es in edgeDict.items():
# 	if len(es) > maxEC:
# 		maxEC = len(es)
# 		maxSeg = s

# print("starting search from ",maxSeg, maxEC)
# longest_path, longestCyclicPath = DFS(maxSeg)

plot_cn_and_multiplicity()

longest_cps = []
longest_path = []
longestCyclicPath = []
for av, cn in scaled_cns.items():
    if cn > 0:
        for v in [av, -av]:
            clp, clcp = DFS(v)
            if len(clp) > len(longest_path):
                longest_path = clp

            if len(clcp) == len(longestCyclicPath) and args.keep_all_LC:
                longest_cps.append(clcp)

            elif len(clcp) > len(longestCyclicPath):
                longestCyclicPath = clcp
                longest_cps = [longestCyclicPath]

total_amp_content = sum([x[2] - x[1] for s, x in id_to_coords.items() if s > 0 and scaled_cns[s] > 0])
print(total_amp_content, " amplified length\n")

print(longest_path, "longest noncyclic")
print(longestCyclicPath, "longest cyclic")


pweights = []
if not args.keep_all_LC:
    paths = [longestCyclicPath, longest_path]

else:
    # remove duplicates
    print("\nremoving duplicate paths\n")
    longest_cps = remove_duplicate_paths(longest_cps)
    longest_cps.append(longest_path)
    paths = longest_cps

plens = []
for p in paths:
    print(p)
    plen = sum([abs(id_to_coords[k][1] - id_to_coords[k][2]) for k in p])
    plens.append(plen)
    print("Path length: " + str(plen))
    cn_remainder_counts = copy.copy(scaled_cns)
    path_amp_content = 0
    seen = set()
    for s in p:
        if abs(s) not in seen:
            path_amp_content += abs(id_to_coords[s][2] - id_to_coords[s][1])
            seen.add(abs(s))

        cn_remainder_counts[abs(s)] -= 1

    print(path_amp_content / total_amp_content, "proportion amplified content explained in path")
    pweights.append(str(path_amp_content / total_amp_content))

    print("remaining unexplained copies")
    for s in sorted(cn_remainder_counts.keys()):
        if cn_remainder_counts[s] > 0:
            print(s, cn_remainder_counts[s])

    print("")

ofname = os.path.basename(args.graph).rsplit("_graph.txt")[0] + "_candidate_cycles.txt"
write_cycles_file(paths, id_to_coords, pweights, scaling_factor, ofname, plens)
