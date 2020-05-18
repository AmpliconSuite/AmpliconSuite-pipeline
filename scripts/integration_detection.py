import sys
import os
import copy
import logging
import pysam as ps
from collections import defaultdict
from intervaltree import Interval, IntervalTree
import argparse
from readclust import *

clustDelta = 500
min_clust_size = 3 #TODO: update this parameter to be based on BAM coverage.
sm_min_clust_size = min_clust_size
padding = 2000 #padding around merged segs to check
tightdiff = 10


def read_excludedRegions(exc_file,ref):
    excIT = defaultdict(IntervalTree)
    with open(exc_file) as infile:
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            fields[1], fields[2] = int(fields[1]), int(fields[2])
            if ref == "GRCh37" and fields[0].startswith("chr"):
                fields[0] = fields[0][3:]

            excIT[fields[0]].add(Interval(fields[1], fields[2]))

    return excIT


def read_graph(graphf):
    gseqs = defaultdict(IntervalTree)
    deList = []
    print("Intervals: ")
    with open(graphf) as infile:
        for line in infile:
            if line.startswith("sequence"):
                fields = line.rstrip().rsplit()
                chrom = fields[1].split(":")[0]
                p1 = int(fields[1].split(":")[1][:-1])
                p2 = int(fields[2].split(":")[1][:-1])
                gseqs[chrom][p1:p2] = float(fields[3])

            if line.startswith("discordant"):
                fields = line.rstrip().rsplit()
                lbp, rbp = fields[1].split("->")
                lchrom, lpd = lbp.rsplit(":")
                rchrom, rpd = rbp.rsplit(":")
                lpos, ldir = int(lpd[:-1]), lpd[-1]
                rpos, rdir = int(rpd[:-1]), rpd[-1]
                rSupp = int(fields[3])
                # isFoldback = (ldir == rdir)

                r1 = dummy_read(lchrom, lpos, ldir == "-")
                r2 = dummy_read(rchrom, rpos, rdir == "-")
                # sr1,sr2 = sorted([r1,r2],key=lambda x: (x.reference_name,x.reference_end))

                curr_clust = pe_read_clust(r1, r2)
                print(str(curr_clust.clust_to_bedpe()))
                curr_clust.size = rSupp
                deList.append(curr_clust)

    return gseqs, deList

#return inSegs,inGraph
def pe_read_in_graph(r1, r2, gseqs, deList):
    chrom1, s1, e1 = r1.reference_name, r1.reference_start, r1.reference_end
    chrom2, s2, e2 = r2.reference_name, r2.reference_start, r2.reference_end
    relSegInts1 = gseqs[chrom1][s1:e1]
    relSegInts2 = gseqs[chrom2][s2:e2]
    inSegs = int(len(relSegInts1) > 0) + int(len(relSegInts2) > 0)
    for gc in deList:
        if gc.rp_has_overlap(r1, r2):
            return inSegs, True

    return inSegs, False


def clust_in_graph(cc, gseqs, deList):
    chrom1, chrom2 = cc.r_IDs
    s1, e1 = cc.centroid
    relSegInts1 = gseqs[chrom1][s1]
    relSegInts2 = gseqs[chrom2][e1]
    inSegs = int(len(relSegInts1) > 0) + int(len(relSegInts2) > 0)
    for gc in deList:
        if gc.clust_has_overlap(cc):
            return inSegs, True

    return inSegs, False


#cn_segs list entry is (chrom,s,e)
def merge_intervals(unsorted_cn_segs):
    msegs = []

    if len(unsorted_cn_segs) < 2:
        return unsorted_cn_segs

    cn_segs = sorted([[x[0],max(0,x[1]-padding), x[2] + padding] for x in unsorted_cn_segs])
    prev_interval = cn_segs[0]
    for i in range(1, len(cn_segs)):
        if prev_interval[0] != cn_segs[i][0] or cn_segs[i][1] - prev_interval[2] > 1:
            msegs.append(prev_interval)
            prev_interval = cn_segs[i]

        else:
            prev_interval[2] = cn_segs[i][2]

    msegs.append(prev_interval)

    return msegs


def clustIsExcludeable(excIT, clust):
    ref1, ref2 = clust.r_IDs[0], clust.r_IDs[0]
    p1, p2 = clust.centroid
    return excIT[ref1][p1] or excIT[ref2][p2]


def readIsExcludeable(excIT,r):
    ref, p = r.reference_name, r.reference_start
    return excIT[ref][p]


#use pysam to extract all of the relevant reads
def get_discordant_reads(alnCollection):
    discordant_read_alns = defaultdict(list)

    for a in alnCollection:

        if not a.is_unmapped and a.is_paired and not a.is_proper_pair and not a.mate_is_unmapped \
                and not a.is_secondary: # and a.mapping_quality >= 5:

            discordant_read_alns[a.query_name].append(a)

    return discordant_read_alns


def sort_filter_discordant_reads(reads,excIT):
    filt_reads = defaultdict(list)

    for k, v in reads.items():

        # ignore reads with multiple alignments (>2). 1 r1, 1 r2
        if len(v) > 2:
            continue

        # one mate in pair falls in the intervals
        elif len(v) == 1:
            #case: handle reads with only one mate mapped
            if v[0].next_reference_id == -1:
                logging.warning("read " + k + " had only one mate in pair mapped")

            #case: one in interval, one out of interval
            dr_r1 = True if v[0].is_read2 else False
            dr_r2 = not dr_r1
            dr = dummy_read(v[0].next_reference_name, v[0].next_reference_start, v[0].mate_is_reverse, v[0].query_name)
            dr.is_read1, dr.is_read2 = dr_r1, dr_r2
            sortedv = (v[0], dr)

        #case: both mates in the intervals
        else:
            if readIsExcludeable(excIT, v[0]) or readIsExcludeable(excIT, v[1]):
                continue

            elif max(v[0].mapping_quality, v[1].mapping_quality) < 15:
                continue

            sortedv = tuple(sorted(v, key=lambda x: (x.reference_name, x.reference_end)))

        refpair = (sortedv[0].reference_name,sortedv[1].reference_name)
        filt_reads[refpair].append(sortedv)

    # sorted_reads = [x for _,x in sorted(zip(r1ends,filt_reads),key=lambda x: x[0])]
    sorted_read_dict = dict()
    for k,l in filt_reads.items():
        sorted_reads = sorted(l,key=lambda x: (x[0].reference_end, x[1].reference_start))
        sorted_read_dict[k] = sorted_reads

    return sorted_read_dict


def cluster_discordant_reads(sr_dict,excIT):
    clusts = defaultdict(list)
    for cp,sr in sr_dict.items():
        curr_clusts = []
        # prev_rIDs = (sr[0][0].reference_name,sr[0][1].reference_name)
        # curr_rIDs = prev_rIDs
        curr_clust = pe_read_clust(sr[0][0],sr[0][1])
        curr_clusts.append(curr_clust)
        for r1,r2 in sr[1:]:
            found = False
            for cc in curr_clusts:
                if r1.reference_end - cc.centroid[0] < 2*clustDelta:
                    if cc.rp_has_overlap(r1,r2):
                        cc.add_pair_to_clust(r1,r2)
                        found = True
                        break

            if not found:
                # rIDs = [r1.reference_name, r2.reference_name]
                curr_clusts.append(pe_read_clust(r1,r2))

        #fencepost at end
        for cc in curr_clusts:
            if cc.size >= min_clust_size:
                if not clustIsExcludeable(excIT,cc):
                    clusts[cp].append(copy.copy(cc))

            elif cc.size >= sm_min_clust_size and cc.total_diff/(cc.size-1) < tightdiff:
                if not clustIsExcludeable(excIT,cc):
                    clusts[cp].append(copy.copy(cc))

    return clusts


def cluster_isLC(cc):
    nonzero_mapq_entries = sum([l.mapping_quality != 0 and r.mapping_quality != 0
                                for l,r in zip(cc.left_reads,cc.right_reads)])

    return nonzero_mapq_entries < 2


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get low-frequency or integration-supporting edges")
    parser.add_argument("--ref", help="Reference genome version.", choices=["hg19", "GRCh37", "GRCh38"], required=True)
    parser.add_argument("--exclude",help="File of breakpoint excludable regions",required=True)
    parser.add_argument("--bam",help="path to bam file",required=True)
    parser.add_argument("--AA_graph",help="path to AA graph",required=True)
    parser.add_argument("-o",help="Output filename prefix")
    args = parser.parse_args()

    base = os.path.basename(args.AA_graph)
    basename = os.path.splitext(base)[0]

    excIT = read_excludedRegions(args.exclude, args.ref)
    if not args.o:
        args.o = basename

    logging.basicConfig(filename=args.o + '_integration.log',
    format = '%(asctime)s %(levelname)-8s %(message)s',
    level = logging.DEBUG,
    datefmt = '%Y-%m-%d %H:%M:%S')

    if not os.path.exists(args.bam):
        sys.stderr.write("BAM file: " + args.bam + " not found\n")
        logging.error("BAM file: " + args.bam + " not found\n")
        sys.exit(1)


    graph_seg_dict, deList = read_graph(args.AA_graph)

    with open(args.o + "_discordant_clusts.txt", 'w') as of1, open(args.o + "_raw_discordant.bedpe", 'w') as of2,\
        open(args.o + "_clust_read_info.txt",'w') as of3:
        of1.write("Sample\tstart_chr\tstart_pos\tend_chr\tend_pos\tNumReads\tinSegs\tinGraph\n")
        of2.write("Sample\tLeftChr\tLeftEnd\tRightChr\tRightStart\tinSegs\tinGraph\n")

        graph_segs = []
        for chrom,it in graph_seg_dict.items():
            for ival in it:
                graph_segs.append([chrom,ival.begin,ival.end])

        logging.debug("graph segs: " + str(graph_segs))
        msegs = merge_intervals(graph_segs)
        logging.debug("msegs: " + str(msegs))


        #iterate over graph segments and get discordant reads.
        logging.info("Extracting reads")
        f_info_dict = {}
        totReads = 0
        allDiscReads = {}
        bamdata = ps.AlignmentFile(args.bam, 'rb')
        relReads = {}
        for seg in msegs:
            chrom, s, e = seg[:3]
            if args.ref == "GRCh37" and chrom.startswith("chr"):
                chrom = chrom[3:]

            alignments = bamdata.fetch(chrom, s, e)
            currNumReads = bamdata.count(chrom, s, e)
            readNamesFreqs = defaultdict(int)
            discordantReads = get_discordant_reads(alignments)
            allDiscReads.update(discordantReads)

        totSings = 0
        for q,l in allDiscReads.items():
            if len(l) == 1:
                totSings+=1

        print(str(totSings) + " single alignment entries")

        logging.info("Sorting and filtering discordant reads")
        sfd_read_dict = sort_filter_discordant_reads(allDiscReads,excIT)
        logging.info("Writing raw read output")
        for cp,rp_l in sfd_read_dict.items():
            for r1,r2 in rp_l:
                if r1.mapping_quality == 0 or r2.mapping_quality == 0:
                    continue
                inSegs, inGraph = pe_read_in_graph(r1,r2, graph_seg_dict, deList)
                of2.write("\t".join([str(x) for x in [basename, r1.reference_name, r1.reference_end, r2.reference_name,
                                                      r2.reference_start, inSegs, inGraph]]) + "\n")

        logging.info("Clustering reads")
        sDR_clusts = cluster_discordant_reads(sfd_read_dict,excIT)

        #filter clusters
        filtered_clusters = defaultdict(list)
        for k,l in sDR_clusts.items():
            for cc in l:
                if not cluster_isLC(cc):
                    filtered_clusters[k].append(cc)

                else:
                    logging.info("removed cluster (Low-Comp): " + str(cc.clust_to_bedpe()))

        #TODO: refine bpoints

        logging.info("Writing read clusts")
        #iterate over clusts
        for k,l in filtered_clusters.items():
            for cc in l:
                #check if a cluster matches something in the graph file.
                inSegs, inGraph = clust_in_graph(cc, graph_seg_dict, deList)
                #write each clust to file
                of1.write("\t".join([str(x) for x in [basename,] + cc.clust_to_bedpe() + [inSegs,inGraph]]) + "\n")
                of3.write(cc.clust_to_string() + "\n")

    logging.info("Finished")
    print("Finished")
    sys.exit()