#!/usr/bin/env python

import os
import sys


# Read the CNVkit .cns files
def convert_cnvkit_cns_to_seeds(cns_file):
    base = os.path.splitext(os.path.basename(cns_file))[0]
    with open(cns_file) as infile, open(base + ".bed", 'w') as outfile:
        head = next(infile).rstrip().rsplit("\t")
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            s, e = int(fields[1]), int(fields[2])
            cn_r = float(fields[4])
            cn = 2 ** (cn_r + 1)
            outline = "\t".join(fields[0:3] + ["CNVkit", str(cn)]) + "\n"
            outfile.write(outline)


if __name__ == '__main__':
    infile = sys.argv[1]
    convert_cnvkit_cns_to_seeds(infile)
