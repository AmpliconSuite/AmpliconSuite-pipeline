from collections import defaultdict

from intervaltree import IntervalTree


# takes list of tuples (chrom, start, end, cn)
def compute_cn_median(cnlist):
    halfn = sum([x[2]-x[1] for x in cnlist])/2
    scns = sorted(cnlist, key=lambda x: x[3])
    rt = 0
    ccn = 0
    for x in scns:
        ccn = x[3]
        rt += (x[2] - x[1])
        if rt > halfn:
            break

    return ccn


# take CNV calls (as bed?) - have to update to not do CNV_GAIN
#input bed file, centromere_dict
#output: path of prefiltered bed file
def prefilter_bed(bedfile, centromere_dict, chr_sizes, cngain, outdir):
    # interval to arm lookup
    region_ivald = defaultdict(IntervalTree)
    for key, value in chr_sizes.items():
        try:
            cent_tup = centromere_dict[key]
            region_ivald[key].addi(0, int(cent_tup[0]), key + "p")
            region_ivald[key].addi(int(cent_tup[1]), int(value), key + "q")

        # handle mitochondrial contig or other things (like viral genomes)
        except KeyError:
            region_ivald[key].addi(0, int(value), key)

    # store cnv calls per arm
    arm2cns = defaultdict(list)
    with open(bedfile) as infile:
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            c, s, e = fields[0], int(fields[1]), int(fields[2]) + 1
            cn = float(fields[-1])
            a = region_ivald[c][(s + e)//2]
            arm2cns[a].append((c, s, e, cn))

    filt_entries = []
    for a in sorted(arm2cns.keys()):
        # compute the median CN of the arm
        init_cns = arm2cns[a]
        med_cn = compute_cn_median(init_cns)
        for x in init_cns:
            ccg = cngain
            if x[2] - x[1] > 10000000:
                ccg *= 1.5

            if x[3] > med_cn + ccg - 2:
                filt_entries.append(x)

    bname = outdir + "/" + bedfile.rsplit("/")[-1].rsplit(".bed")[0] + "_pre_filtered.bed"
    with open(bname, 'w') as outfile:
        outfile.write("\t".join([str(x) for x in filt_entries]) + "\n")

    return bname
