from collections import defaultdict
import subprocess
import sys

# create the set of autosomal chromosome names in various builds.
# should be updated if a reference is added to the data repo with more than 22 autosomes, but not necessary to do so
chrom_range = [str(x) for x in range(1, 23)]
chrom_range.extend(["chr" + x for x in chrom_range])
chrom_range = set(chrom_range)


# get a subset of the chromosome names/lengths from a .fai file.
def get_ref_seq_lens(ref_genome_size_file):
    chr_sizes = {}
    try:
        with open(ref_genome_size_file) as infile:
            for line in infile:
                fields = line.rstrip().rsplit()
                if fields[0] in chrom_range:
                    chr_sizes[fields[0]] = int(fields[1])

    except IOError:
        pass

    return chr_sizes


# read bam header and store info
def get_bam_header(bamf):
    cmd = 'samtools view -H ' + bamf
    return subprocess.check_output(cmd, shell=True).decode("utf-8")


# extract sequence lengths and ids
def extract_seq_info(bam_header):
    bamSeqLenD = defaultdict(int)
    linelist = bam_header.rsplit("\n")
    for line in (x for x in linelist if x.startswith("@SQ")):
        fields = line.rstrip().rsplit()[1:]
        ld = {i.rsplit(":")[0]: i.rsplit(":")[1] for i in fields}
        bamSeqLenD[ld["SN"]] = int(ld["LN"])

    return bamSeqLenD


# check if bam matches to a reference genome in terms of length and sequence name
# returns false if the same chromosome has different length in bam vs. reference
# returns false if no chromosome names are shared between bam/reference
# returns true if no shared chromosomes have different lengths and at least one chromosome is present.
def match_ref(bamSeqLenD, ref_len_d):
    overlaps = False
    for chrom, len in ref_len_d.items():
        if bamSeqLenD[chrom] > 0 and len != bamSeqLenD[chrom]:
            return False

        elif len == bamSeqLenD[chrom]:
            overlaps = True

    return overlaps


# check if the BAM reference matches to sequence names & lengths in a dictionary of .fai files
# returns the name of the reference genome the BAM matches to, or prints error and returns None.
def check_ref(bamf, ref_to_fai_dict):
    bam_header = get_bam_header(bamf)
    bamSeqLenD = extract_seq_info(bam_header)
    for refName, fai_path in ref_to_fai_dict.items():
        ref_len_d = get_ref_seq_lens(fai_path)
        matched = match_ref(bamSeqLenD, ref_len_d)
        if matched:
            print("Matched " + bamf + " to reference genome " + refName)
            return refName

    sys.stderr.write("ERROR: Could not match BAM to a known AA reference genome!\n")
    sys.stderr.write("This may happen if 1) The value provided to optional argument '--ref' does not match the "
                     "reference the BAM is aligned to, or 2) The corresponding AA data repo folder for this reference "
                     "is not present, or 3) The BAM uses an obscure chromosome naming convention. Consider inspecting "
                     "the header of the BAM file and the AA data repo directory.\n")

    return None