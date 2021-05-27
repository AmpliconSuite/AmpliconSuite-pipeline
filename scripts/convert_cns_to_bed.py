#!/usr/bin/env python

import os
import sys

# takes one argument as input, a .cns file from CNVKit, rewrites it as a bed after converting the CN data.

# Read the CNVkit .cns files
def convert_cnvkit_cns_to_seeds(cns_file, base):
    ofname = base + "_uncorr_CN.bed"
    with open(cns_file) as infile, open(ofname, 'w') as outfile:
        head = next(infile).rstrip().rsplit("\t")
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            # s, e = int(fields[1]), int(fields[2])
            cn_r = float(fields[4])
            cn = 2 ** (cn_r + 1)
            outline = "\t".join(fields[0:3] + ["CNVkit", str(cn)]) + "\n"
            outfile.write(outline)

    return ofname


def estimate_ploidy_correction(uncorr, base):
    with open(uncorr) as infile:
        tl = 0
        tw = 0.0
        allw = []
        for line in infile:
            fields = line.rstrip().rsplit()
            cname = fields[0].lstrip('chr')
            if cname != 'X' and cname != 'Y' and cname != 'M':
                l = int(fields[2]) - int(fields[1])
                w = float(fields[4]) * l
                allw.append(float(fields[4]))

                tl += l
                tw += w

    # print((tw / tl) / 2.0)
    r = (tw / tl) / 2.0

    with open(uncorr) as infile, open(base + "_ESTIMATED_PLOIDY_CORRECTED_CN.bed", 'w') as outfile:
        for line in infile:
            fields = line.rstrip().rsplit()
            outfile.write("\t".join(fields[:4]))
            outfile.write("\t" + str(float(fields[4]) * r) + "\n")


if __name__ == '__main__':
    infile = sys.argv[1]
    base = os.path.splitext(os.path.basename(infile))[0]
    uncorr = convert_cnvkit_cns_to_seeds(infile, base)
    estimate_ploidy_correction(uncorr, base)
