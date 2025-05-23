from collections import defaultdict
import logging
import subprocess
import sys

# create the set of autosomal chromosome names in various builds.
# should be updated if a reference is added to the data repo with more than 22 autosomes, but not necessary to do so
chrom_range = [str(x) for x in range(1, 23)]
chrom_range.extend(["chr" + x for x in chrom_range])
chrom_range.append("hpv16ref_1")  # use one representative entry from the viral genome collection to catch a viral ref.
chrom_range = set(chrom_range)


def get_ref_fname(aa_dr_path, rname):
    with open(aa_dr_path + "/" + rname + "/file_list.txt") as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            if fields[0] == "fa_file":
                return fields[1]

    logging.error("ERROR: AA data repo 'file_list.txt' not found!\n")
    return None


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
def get_bam_header(bamf, samtools):
    cmd = samtools + ' view -H ' + bamf
    return subprocess.check_output(cmd, shell=True).decode("utf-8")


# extract sequence lengths and ids
def extract_seq_info(bam_header):
    bamSeqLenD = defaultdict(int)
    linelist = bam_header.rsplit("\n")
    altFound = False
    for line in (x for x in linelist if x.startswith("@SQ")):
        fields = line.rstrip().rsplit()[1:]
        ld = {i.rsplit(":")[0]: i.rsplit(":")[1] for i in fields if ":" in i}
        bamSeqLenD[ld["SN"]] = int(ld["LN"])
        if ld["SN"].endswith("_alt"):
            altFound = True

    if altFound:
        logging.warning("WARNING: Alt contigs detected in bam header. If BAM alignment was not alt-aware, AA results may be incorrect!\n")

    return bamSeqLenD


""" 
check if bam matches to a reference genome in terms of length and sequence name
returns false if the same chromosome has different length in bam vs. reference
returns false if no chromosome names are shared between bam/reference
returns true if no shared chromosomes have different lengths and at least one chromosome is present.
"""
def match_ref(bamSeqLenD, ref_len_d):
    overlaps = 0
    for chrom, len in ref_len_d.items():
        if bamSeqLenD[chrom] > 0 and len != bamSeqLenD[chrom]:
            return False

        elif len == bamSeqLenD[chrom]:
            overlaps+=1

    return overlaps


# check properly paired rate on bam file
def check_properly_paired(bamf, samtools):
    cmd = samtools + " flagstat {} | grep 'properly paired'".format(bamf)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    exit_code = process.returncode

    # Convert stdout to a string (in Python 3 it is already a string, but in Python 2 it is bytes)
    t = stdout.decode("utf-8") if isinstance(stdout, bytes) else stdout

    logging.info(bamf + ": " + t.rstrip() + "\n")
    if int(exit_code) != 0:
        em = stderr.decode("utf-8") if isinstance(stderr, bytes) else stderr
        logging.error(em)
        logging.error("Samtools flagstat returned a non-zero exit code or reported no information on properly paired reads! "
                      "This indicates a significant problem with your bam file.")
        sys.exit(1)

    ppp = float(t.rsplit("(")[-1].rsplit("%")[0])
    if t.startswith("0 + 0"):
        logging.error("ERROR: UNSUITABLE BAM FILE! No properly-paired reads were found by samtools. "
                         "AmpliconSuite-pipeline requires paired-end sequencing data. AA requires paired-end WGS. If "
                      "this was PE WGS, please confirm your alignment steps were done appropriately\n\n")
        sys.exit(1)

    elif ppp < 95:
        logging.warning("WARNING: BAM FILE PROPERLY PAIRED RATE IS BELOW 95%.\nQuality of data may be insufficient for AA "
              "analysis. Poorly controlled insert size distribution during library prep can cause high fractions of read"
              " pairs marked as discordant during alignment. Artifactual short SVs and long runtimes may occur!"
              "\n")

    return ppp


# check if the BAM reference matches to sequence names & lengths in a dictionary of .fai files
# returns the name of the reference genome the BAM matches to, or prints error and returns None.
def check_ref(bamf, ref_to_fai_dict, samtools):
    bam_header = get_bam_header(bamf, samtools)
    bamSeqLenD = extract_seq_info(bam_header)
    bestref = None
    bestrefhits = 0
    for refName, fai_path in ref_to_fai_dict.items():
        ref_len_d = get_ref_seq_lens(fai_path)
        matched = match_ref(bamSeqLenD, ref_len_d)
        if matched:
            if matched > bestrefhits:
                bestref = refName
                bestrefhits = matched

            elif bestref and matched == bestrefhits and "_viral" in bestref and "_viral" not in refName:
                bestref = refName
                bestrefhits = matched

    if bestref:
        logging.info("Matched " + bamf + " to reference genome " + bestref)
        return bestref

    em1 = "ERROR: Could not match BAM to a known AA reference genome!\n"
    em2 = """This may happen if 1) The value provided to optional argument '--ref' does not match the 
          reference the BAM is aligned to, or 2) The corresponding AA data repo folder for this reference 
          is not present, or 3) The BAM uses a different chromosome naming convention (e.g. accession 
          numbers instead of chromosome names). Consider inspecting the header of the BAM file and the AA data repo 
          directory. Additional data repo reference genomes can be acquired by running with --download_repo [ref]\n"""

    logging.error(em1)
    logging.error(em2)
    # sys.stderr.write(em1)
    # sys.stderr.write(em2)

    return None
