import sys
import os
import copy
import logging
import pysam as ps
from collections import defaultdict
from intervaltree import Interval, IntervalTree
import argparse

clustDelta = 500
min_clust_size = 25
fb_dist_cut = 25000

#for bfb-calling
min_score_for_bfb = 0.25
min_fb_reads_for_bfb = 25

class dummy_read(object):
    def __init__(self,mrn,mstart,mir):
        self.reference_start = mstart
        self.reference_id = mrn
        self.is_reverse = mir
        self.qstart,self.qend = -1,-1
        self.reference_end = mstart
        self.template_length = -1
        self.is_read1 = False
        self.is_read2 = True


#assume everything from same chromosome
class read_clust(object):
    def __init__(self,rname,r1,r2,rIDs,is_foldback=False):
        self.clust_ID = None
        self.left_reads, self.right_reads, self.r_IDs = [],[],[]
        self.size = 0
        self.rIDs = rIDs
        self.centroid = (0,0)
        self.is_foldback = is_foldback
        self.add_pair_to_clust(rname,r1,r2)

    def add_pair_to_clust(self,rname,r1,r2):
        # self.clust_ID = rname
        self.left_reads.append(r1)
        self.right_reads.append(r2)
        self.r_IDs.append(rname)
        self.size+=1
        self.update_centroid()
        self.clust_ID = str(self.centroid)

    def update_centroid(self):
        currL,currR = self.centroid
        prevSS = len(self.left_reads)-1
        wL = prevSS*currL
        wR = prevSS*currR
        meanL = (wL + self.left_reads[-1].reference_end)/len(self.left_reads)
        meanR = (wR + self.right_reads[-1].reference_start)/len(self.right_reads)
        self.centroid = (meanL,meanR)

    def has_overlap(self,r1,r2):
        if r1.reference_id != self.rIDs[0] or r2.reference_id != self.rIDs[1]:
            return False

        return abs(r1.reference_end - self.centroid[0]) < clustDelta and \
               abs(r2.reference_start - self.centroid[1]) < clustDelta

    def clust_to_bedpe(self):
        a = self.left_reads[-1]
        b = self.right_reads[-1]
        return [str(a.reference_id),str(self.centroid[0]),str(b.reference_id),str(self.centroid[1]),
                          str(int(self.is_foldback)),str(self.size)]

    # def clust_to_string(self):
    #     s = self.clust_ID + " | #read_pairs: " + str(self.size) + "\n"
    #
    #     for v in zip(self.left_reads,self.right_reads):
    #         for a in v:
    #             readno = ""
    #             if a.is_read1:
    #                 readno = "Read 1:  "
    #             elif a.is_read2:
    #                 readno = "Read 2:  "
    #
    #             adir = ""
    #             if a.is_reverse:
    #                 a_dir = "-"
    #             else:
    #                 a_dir = "+"
    #
    #             sstart = readno + ": " + adir + " "
    #             s+=sstart
    #             s+=" ".join([str(x) for x in [a.qstart,a.qend,a.reference_id,a.reference_start,a.reference_end,
    #                                          a.is_reverse,a.template_length]])
    #             s+="\n"
    #
    #         s+="\n"
    #
    #     s+="\n"
    #     return s


def read_excludedRegions(exc_file,ref):
    excIT = defaultdict(IntervalTree)
    with open(exc_file) as infile:
        for line in infile:
            line = line.rstrip()
            if not line:
                continue

            fields = line.rsplit("\t")
            fields[1],fields[2] = int(fields[1]),int(fields[2])
            if ref == "GRCh37" and fields[0].startswith("chr"):
                fields[0] = fields[0][3:]

            excIT[fields[0]].add(Interval(fields[1],fields[2]))

    return excIT


def parse_bfb_file(fname):
    with open(fname) as infile:
        _ = next(infile)
        outputline = next(infile)
        vect = [int(x) for x in outputline.rsplit(" (")[0][1:-1].rsplit(", ")]

        return vect


def parse_cnv_file(fname):
    tuplist = []
    with open(fname) as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            fields[1] = int(fields[1])
            fields[2] = int(fields[2])
            fields[3] = int(fields[3])
            tuplist.append(tuple(fields))

    return tuplist


def filter_and_merge_intervals(cn_vect, cn_segs, is_parm):
    #filter
    filt_segs = []
    fm_segs = []

    if max(cn_vect) > 1 and len(set(cn_vect)) > 1:
        for c, x in zip(cn_vect, cn_segs):
            if c >= 1:
                filt_segs.append(x)

        if is_parm:
            filt_segs = filt_segs[::-1]

        # merge
        # assume more than three
        curr_interval = list(filt_segs[0])
        for i in range(1, len(filt_segs)):
            if filt_segs[i][1] - curr_interval[2] <= 1:
                curr_interval[2] = filt_segs[i][2]

            else:
                fm_segs.append(curr_interval)
                curr_interval = list(filt_segs[i])

        fm_segs.append(curr_interval)

    return fm_segs


def get_discordant_reads(alnCollection):
    discordant_alns = {}
    for a in alnCollection:
        if not a.is_unmapped and a.is_paired and not a.is_proper_pair and not a.mate_is_unmapped \
                and not a.is_secondary and a.mapping_quality >= 5:

            if a.query_name not in discordant_alns:
                discordant_alns[a.query_name] = []

            discordant_alns[a.query_name].append(a)

    return discordant_alns


def isExcludeable(excIT,clust):
    ref1,ref2 = clust.rIDs[0],clust.rIDs[0]
    p1,p2 = clust.centroid
    return excIT[ref1][p1] or excIT[ref2][p2]


def sort_filter_discordant_reads(reads):
    filt_reads, r1ends = [],[]

    for k,v in reads.items():
        if len(v) > 2:
            continue
        elif len(v) == 1:
            if v[0].next_reference_id == -1:
                print("issue")
            dr = dummy_read(v[0].next_reference_id,v[0].next_reference_start, not v[0].is_reverse)
            sortedv = [v[0],dr]

        else:
            sortedv = sorted(v,key=lambda x: x.reference_end)

        filt_reads.append(sortedv)
        r1ends.append(sortedv[0].reference_end)

    sorted_reads = [x for _,x in sorted(zip(r1ends,filt_reads),key=lambda x: x[0])]
    sDR,sFB = [],[]
    for l,r in sorted_reads:
        if (l.is_reverse == r.is_reverse):
            if l.reference_id == r.reference_id and abs(r.reference_start - l.reference_end) < fb_dist_cut:
                sFB.append((l.query_name,l,r))
            else:
                sDR.append((l.query_name,l,r))
        else:
            sDR.append((l.query_name,l,r))

    return sDR,sFB


def cluster_discordant_reads(sr,excIT):
    clusts = []
    curr_clusts = []
    if sr:
        rIDs = [sr[0][1].reference_id,sr[0][2].reference_id]
        curr_clust = read_clust(sr[0][0],sr[0][1],sr[0][2],rIDs)
        curr_clusts.append(curr_clust)
        for rname,r1,r2 in sr[1:]:
            popCount = 0
            found = False
            for cc in curr_clusts:
                if r1.reference_end - cc.centroid[0] < 2*clustDelta:
                    if cc.has_overlap(r1,r2):
                        cc.add_pair_to_clust(rname,r1,r2)
                        found = True
                        break

                else:
                    popCount+=1

            if not found:
                rIDs = [r1.reference_id, r2.reference_id]
                curr_clusts.append(read_clust(rname,r1,r2,rIDs))

            for i in range(popCount):
                cc = curr_clusts[i]
                if cc.size >= min_clust_size:
                    if not isExcludeable(excIT,cc):
                        clusts.append(copy.copy(curr_clusts[i]))

            curr_clusts = curr_clusts[popCount:]

        #fencepost at end
        for cc in curr_clusts:
            if cc.size >= min_clust_size:
                if not isExcludeable(excIT,cc):
                    clusts.append(copy.copy(cc))

    return clusts


#Compute f from the edges in the AA graph alone
def compute_f_from_AA_graph(graphf,excIT):
    fbCount,nonFbCount = 0,0
    with open(graphf) as infile:
        for line in infile:
            if line.startswith("discordant"):
                fields = line.rstrip().rsplit()
                lbp,rbp = fields[1].split("->")
                lchrom,lpd = lbp.rsplit(":")
                rchrom,rpd = rbp.rsplit(":")
                lpos,ldir = int(lpd[:-1]),lpd[-1]
                rpos,rdir = int(rpd[:-1]),rpd[-1]
                rSupp = int(fields[3])
                if ldir == rdir:
                    #check if goes to forbidden
                    if excIT[lchrom][lpos] or excIT[rchrom][rpos]:
                        continue

                    elif lchrom == rchrom and abs(rpos - lpos) < fb_dist_cut:
                        fbCount+=rSupp

                    else:
                        nonFbCount+=rSupp

                else:
                    nonFbCount+=rSupp

    return fbCount,nonFbCount


if __name__ == '__main__':
    # data_table_fname = sys.argv[1]
    # data_dir = sys.argv[2]
    # data_dir = "/nucleus/projects/namphuon/cell_line/chromothripsis/"
    # bam_dir = "/nucleus/pedigree/projects/extrachromosome/data/turner2017/bwa_bam/"
    # data_table_fname = "BFB_calls.xlsx"

    parser = argparse.ArgumentParser(description="Compute f for a single sample")
    parser.add_argument("--bam",help="path to bam file")
    parser.add_argument("--exclude",help="File of breakpoint excludable regions",required=True)
    parser.add_argument("-o",help="Output filename prefix")
    parser.add_argument("--ref", help="Reference genome version.", choices=["hg19", "GRCh37", "GRCh38"], default="GRCh37")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--cnv_dir",help="path to samples directory with CNV & BFB results")
    group.add_argument("--AA_graph_list",help="path to list of AA graphs")
    args = parser.parse_args()

    if args.cnv_dir:
        if not args.bam or not args.exclude or not args.ref:
            print("--bam, --ref, and --exclude required with --cnv_dir")
            sys.exit(1)

        basename = args.cnv_dir.rsplit("/")[-1] if not args.cnv_dir.endswith("/") else args.cnv_dir[:-1].rsplit("/")[-1]
        # read excludable intervals
    else:
        t = args.AA_graph_list
        #basename = t.rsplit("/")[-1] if not t.endswith("/") else t[:-1].rsplit("/")[-1]
        base = os.path.basename(t)
        basename = os.path.splitext(base)[0]

    excIT = read_excludedRegions(args.exclude, args.ref)
    if not args.o:
        args.o = basename

    logging.basicConfig(filename=args.o + '.log',
    format = '%(asctime)s %(levelname)-8s %(message)s',
    level = logging.DEBUG,
    datefmt = '%Y-%m-%d %H:%M:%S')



    if args.bam and not os.path.exists(args.bam):
        sys.stderr.write("BAM file: " + args.bam + " not found\n")
        logging.error("BAM file: " + args.bam + " not found\n")
        sys.exit(1)

    # If BFB strings provided
    if args.cnv_dir:
        logging.info("Reading BFB strings and CN data")
        #list the directory
        flist = os.listdir(args.cnv_dir)
        cnvFiles = {}
        bfbFiles = {}
        #pair files
        for f in flist:
            if f.endswith(".cnv"):
                arm,pos = f.rsplit(".cnv")[0].rsplit('segments.')[1].rsplit('.')
                cnvFiles[(arm,pos)] = args.cnv_dir + "/" + f

            elif f.startswith("bfb."):
                arm,pos = f.rsplit("bfb.")[1].rsplit(".")
                bfbFiles[(arm,pos)] = args.cnv_dir + "/" + f

        id_to_regions = {}

        with open(args.o + "_clustf.txt", 'w') as of1, open(args.o + "_rawclust.bedpe", 'w') as of2:
            of1.write("Sample\tChr\tArm\tClust-f\tNumFB\tTotDiscordantInClust\n")
            of2.write("Sample\tChr\tArm\tLeftChr\tLeftPos\tRightChr\tRightPos\tIsFoldback\tnReads\n")
            for chrom, arm in cnvFiles:
                chromarm = "." + chrom + "." + arm
                bfbVect = parse_bfb_file(bfbFiles[(chrom, arm)])
                cn_data = parse_cnv_file(cnvFiles[(chrom, arm)])
                id_to_regions[(basename, chrom, arm)] = [bfbVect, cn_data]
                fm_segs = filter_and_merge_intervals(bfbVect, cn_data, arm == "p")
                if not fm_segs:
                    of1.write("\t".join([basename, chrom, arm, "invalid"]) + "\n")
                    continue

                logging.info(basename + " " + chrom + " " + arm + ": Extracting reads and computing f...")
                # now use pysam to extract all of the relevant reads
                f_info_dict = {}
                samp_clust_vect = [[]] * len(id_to_regions)
                totReads = 0
                allDiscReads = {}
                bamdata = ps.AlignmentFile(args.bam, 'rb')
                relReads = {}
                nFB = 0
                total_bfb_amp_len = 0
                for seg in fm_segs:
                    chrom, s, e = seg[:3]
                    if args.ref == "GRCh37" and chrom.startswith("chr"):
                        chrom = chrom[3:]

                    total_bfb_amp_len += (e - s)
                    alignments = bamdata.fetch(chrom, s, e)
                    currNumReads = bamdata.count(chrom, s, e)
                    readNamesFreqs = defaultdict(int)
                    discordantReads = get_discordant_reads(alignments)
                    allDiscReads.update(discordantReads)

                # print(len(allDiscReads), "discordant reads")
                sDR, sFB = sort_filter_discordant_reads(allDiscReads)

                logging.info("Clustering reads")
                # print(len(sDR),len(sFB),"sizes of sDR and sFB")
                sDR_clusts = cluster_discordant_reads(sDR, excIT)
                sFB_clusts = cluster_discordant_reads(sFB, excIT)

                logging.info("Writing output\n")
                totR = 0
                totFB = 0
                for cc in sDR_clusts:
                    of2.write("\t".join([basename, chrom, arm] + cc.clust_to_bedpe()) + "\n")
                    totR += cc.size

                for cc in sFB_clusts:
                    of2.write("\t".join([basename, chrom, arm] + cc.clust_to_bedpe()) + "\n")
                    totR += cc.size
                    totFB += cc.size

                totR = max(totR, 1)
                clustF = (totFB / totR)
                numnonfb = (totR - 1) - totFB
                strVals = [str(x) for x in [clustF, totFB, numnonfb]]
                outstring = "\t".join([basename, chrom, arm] + strVals) + "\n"
                of1.write(outstring)


    elif args.AA_graph_list:
        samp_list = []
        samp_with_bfb = set()
        with open(args.AA_graph_list) as infile, open(args.o + "_graph_f.txt", 'w') as outfile:
            outfile.write("Sample\tgraph-f\tfoldback_reads\tnon_foldback_reads\n")
            for line in infile:
                try:
                    a,f = line.rstrip().split()
                    totFB, totNonFB = compute_f_from_AA_graph(f,excIT)
                    gF = float(totFB) / max(1,(totFB + totNonFB))

                    samp = a.rsplit("_amplicon")[0]
                    if gF >= min_score_for_bfb and totFB >= min_fb_reads_for_bfb:
                        samp_with_bfb.add(samp)

                    if samp not in samp_list:
                        samp_list.append(samp)

                except ValueError:
                    a = line.rstrip()
                    samp = a.rsplit("_amplicon")[0]
                    if samp not in samp_list:
                        samp_list.append(samp)

                    totFB = 0
                    totNonFB = 0
                    gF = 0

                outfile.write("\t".join([a,str(gF),str(totFB),str(totNonFB)]) + "\n")

        #write bfb calls
        with open(args.o + "_bfb_calls.txt", 'w') as outfile:
            for samp in samp_list:
                isBFB = str(samp in samp_with_bfb)
                outfile.write(samp + "\t" + isBFB + "\n")

    sys.exit()






