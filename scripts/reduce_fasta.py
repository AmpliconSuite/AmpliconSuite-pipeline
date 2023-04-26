#!/usr/bin/env python

import argparse
from itertools import groupby
import os
import sys

# reduce fasta to only entries in given list or .fai.


# read list of chroms from .fai or file of one chrom name per line.
# strips ">"
def getRelChrs(CHROMS):
    chrSet = set()
    with open(CHROMS) as chromFile:
        for line in chromFile:
            line = line.rstrip()
            if line:
                chrName = line.rsplit()[0].lstrip(">")
                chrSet.add(chrName)

    chrList = sorted(list(chrSet))
    return chrList


def fasta_reader(fasta_file, chroms_to_get):
    fasta_dict = {}
    print("Reading FASTA: {}".format(fasta_file))
    with open(fasta_file) as infile:
        faiter = (x[1] for x in groupby(infile, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            seq_name = next(header)[1:].rstrip().rsplit()[0]
            if (seq_name in chroms_to_get):
                # join all sequence lines to one.
                seq = "".join(s.strip() for s in next(faiter))
                fasta_dict[seq_name] = seq

    return fasta_dict


if __name__ == '__main__':
    # Parses the command line arguments
    parser = argparse.ArgumentParser(description="Reduce FASTA to entries in list")
    parser.add_argument("-r", "--ref", help="reference fasta",required=True)
    parser.add_argument("-c", "--chrom", help="reference fasta chrom list file. assumes one entry per line",required=True)
    parser.add_argument("-o", "--outname", help="output file directory")

    args = parser.parse_args()
    if not args.outname:
        args.outname = ""

    chrList = getRelChrs(args.chrom)
    seqD = fasta_reader(args.ref, chrList)
    base = os.path.basename(args.ref)
    refGOutName = os.path.splitext(base)[0] + "_reduced" + "".join(os.path.splitext(base)[1:])
    print("Writing stripped FASTA\n")
    with open(args.outname + refGOutName, 'w') as outfile:
        for i in chrList:
            outfile.write(">" + i + "\n")
            outfile.write(seqD[i] + "\n")

    sys.exit()
