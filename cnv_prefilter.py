from collections import defaultdict
import logging
import os

from intervaltree import IntervalTree


def merge_intervals(usort_intd, cn_cut=4.5, tol=1, require_same_cn=False):
    merged_intd = defaultdict(IntervalTree)
    for chrom, usort_ints in usort_intd.items():
        # sort ints
        sort_ints = sorted([x for x in usort_ints if x[2] > cn_cut])
        if not sort_ints:
            continue

        # merge sorted ints
        mi = [sort_ints[0]]
        for ival in sort_ints[1:]:
            pass_cn_check = True
            if require_same_cn and not ival[2] == mi[-1][2]:
                pass_cn_check = False
                
            if ival[0] <= mi[-1][1] + tol and pass_cn_check:
                ui = (mi[-1][0], max(ival[1], mi[-1][1]), mi[-1][2])
                mi[-1] = ui

            else:
                mi.append(ival)

        for x in mi:
            merged_intd[chrom].addi(x[0], x[1], x[2])

    return merged_intd


# create an interval list (chrom, start, end, CN) from a dict of interval trees.
def ivald_to_ilist(ivald):
    ivals = []
    for chrom, ivalt in ivald.items():
        for ival in ivalt:
            ivals.append((chrom, ival.begin, ival.end, ival.data))

    return ivals


# takes list of tuples (chrom, start, end, cn)
def compute_cn_median(cnlist, armlen):
    cnsum = sum([x[2]-x[1] for x in cnlist])
    if cnsum < 0.5 * armlen:
        return 2.0

    halfn = cnsum/2.0
    scns = sorted(cnlist, key=lambda x: x[3])
    rt = 0
    ccn = 0
    for x in scns:
        ccn = x[3]
        rt += (x[2] - x[1])
        if rt >= halfn:
            break

    return ccn


def read_bed(ifname, keepdat=False):
    beddict = defaultdict(IntervalTree)
    with open(ifname) as infile:
        for line in infile:
            line = line.rstrip()
            if line:
                fields = line.rsplit()
                s, e = int(fields[1]), int(fields[2])
                if e - s == 0:
                    logging.warning("Size 0 interval found. Skipping: " + line)
                    continue

                if keepdat:
                    beddict[fields[0]].addi(s, e, tuple(fields[3:]))
                else:
                    beddict[fields[0]].addi(s, e)

    return beddict


# read regions to split on/filter into dictionary of interval trees, where keys are chromosomes
def read_gain_regions(ref):
    AA_DATA_REPO = os.environ["AA_DATA_REPO"] + "/" + ref + "/"
    fdict = {}
    with open(AA_DATA_REPO + "file_list.txt") as infile:
        for line in infile:
            line = line.rstrip()
            if line:
                fields = line.rsplit()
                fdict[fields[0]] = fields[1]

    grf = AA_DATA_REPO + fdict["conserved_regions_filename"]
    gain_regions = read_bed(grf)

    return gain_regions


def get_continuous_high_regions(bedfile, cngain):
    raw_input = defaultdict(list)
    with open(bedfile) as infile:
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            c, s, e = fields[0], int(fields[1]), int(fields[2]) + 1
            cn = float(fields[-1])
            raw_input[c].append((s,e,cn))

    return merge_intervals(raw_input, cn_cut=cngain, tol=300000)


# take CNV calls (as bed?) - have to update to not do CNV_GAIN
#input bed file, centromere_dict
#output: path of prefiltered bed file
def prefilter_bed(bedfile, ref, centromere_dict, chr_sizes, cngain, outdir):
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
    arm2lens = defaultdict(int)
    with open(bedfile) as infile:
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            c, s, e = fields[0], int(fields[1]), int(fields[2]) + 1
            if c == "hs37d5":
                continue

            cn = float(fields[-1])
            a = region_ivald[c][(s + e)//2]
            if not a:
                a = region_ivald[c][s:e]

            if a:
                carm_interval = a.pop()
                carm = carm_interval.data
                arm2cns[carm].append((c, s, e, cn))
                arm2lens[carm] = carm_interval.end - carm_interval.begin

            else:
                arm2cns["other"].append((c, s, e, cn))
                # print("Warning: could not match " + c + ":" + str(s) + "-" + str(e) + " to a known chromosome arm!")

    continuous_high_region_ivald = get_continuous_high_regions(bedfile, cngain)
    cn_filt_entries = []
    for a in sorted(arm2cns.keys()):
        # compute the median CN of the arm
        init_cns = arm2cns[a]
        med_cn = compute_cn_median(init_cns, arm2lens[a])
        for x in init_cns:
            ccg = cngain
            continuous_high_hits = continuous_high_region_ivald[x[0]][x[1]:x[2]]
            if continuous_high_hits:
                for y in continuous_high_hits:
                    if y.end - y.begin > 10000000:
                        ccg *= 1.5
                        break

            if x[3] > med_cn + ccg - 2:
                cn_filt_entries.append(x)
            elif ref == "GRCh38_viral" and not x[0].startswith("chr") and x[3] > 0.1:
                cn_filt_entries.append(x)

    gain_regions = read_gain_regions(ref)
    # now remove regions based on filter regions
    filt_ivald = defaultdict(IntervalTree)
    for x in cn_filt_entries:
        cit = IntervalTree()
        cit.addi(x[1], x[2])
        bi = gain_regions[x[0]]
        for y in bi:
            cit.slice(y.begin)
            cit.slice(y.end)

        for p in sorted(cit):
            filt_ivald[x[0]].addi(p[0], p[1], x[3])

    merged_filt_ivald = merge_intervals(filt_ivald, cn_cut=cngain, require_same_cn=True)
    final_filt_entries = ivald_to_ilist(merged_filt_ivald)
    bname = outdir + "/" + bedfile.rsplit("/")[-1].rsplit(".bed")[0] + "_pre_filtered.bed"
    with open(bname, 'w') as outfile:
        for entry in final_filt_entries:
            outfile.write("\t".join([str(x) for x in entry]) + "\n")

    return bname
