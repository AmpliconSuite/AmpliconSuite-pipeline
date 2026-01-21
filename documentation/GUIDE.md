## Primer on running AmpliconArchitect

### Contents
 - [What is AA?](#what-is-aa)
 - [Preparing the inputs](#preparing-the-inputs)
 - [Selecting seed regions](#filtering-and-selecting-seeds-from-cnv-data)
 - [Worked example](#worked-example)
 - [FAQ](#faq)

#

### What is AA?
AmpliconArchitect (AA) is a tool designed to study **focal amplifications** in the cancer genome. AA can help identify common sources of focal amplifications, such as **circular extrachromosomal DNA (ecDNA)**, and **breakage-fusion-bridge cycles (BFBs)**.
 
AA takes as input WGS data and a bed file of regions to examine (called "seed regions"), and outputs a **breakpoint graph** and a **cycles file**. The **breakpoint graph**
encodes the identity of regions in the genome, their copy numbers (as measured by AA), and the breakpoint junctions between
these elements (as measured by AA). The **cycles file** outputs both non-cyclic paths and cyclic paths which are decompositions
of the breakpoint graph into paths which best explain the CN of the segments.


An open-access publication detailing AA's methodology is available in [Nature Communications](https://www.nature.com/articles/s41467-018-08200-y). Please cite AA if you use it in your work. 

### How does AA work?
Given the seed regions, AA performs multiple cycles of exploration for additional SVs joining the graph regions to other locations in the genome and recruits those genome intervals to the graph, ultimately forming what we termed an “amplicon.” The amplicon graph contains three general types of SV edges, “concordant edges” that connect directly adjacent segments of the genome, “discordant edges” that connect non-adjacent pieces of the genome [colored by orientation](https://github.com/virajbdeshpande/AmpliconArchitect?tab=readme-ov-file#4-the-sv-view-out_ampliconidpngpdf), and “source edges” that exit the graph segments to neighboring locations not included in the graph, or which represent SVs with one end located inside the graph and the other end being an unknown location. AA applies a balanced-flow constraint to correct copy-numbers of genomic segments and SV edges, which it solves using convex optimization. AA then explores the genome graph to decompose it into paths and cycles, constrained by the copy-number available. For focal amplifications of simple structure, individual decompositions may capture the full structure of the focal amplification, but in complex cases they typically represent substructures of a larger or more heterogeneous amplicon.

#

### What is AmpliconSuite-pipeline?
AmpliconSuite-pipeline is a workflow that runs AmpliconArchitect, as well as upstream steps (alignment, seed region identification) and downstream steps (AmpliconClassifier, packaging outputs for AmpliconRepository).

### Preparing the inputs
![AA workflow](../images/AA_example.png)

AA takes as input a WGS BAM file (paired-end WGS), and a user-created BED file of seed regions as inputs. Here we will discuss some of 
the best practices for generating these files.

AA uses external CNV calls to determine which regions it should examine - thes are called **CNV seeds**. 
However, AA independently calls copy number inside the seed regions it is tasked with, and thus after selecting the regions, 
**those CN calls are not propogated into AA's own estimations.**

**To help standardize the process of running AA, we have created a wrapper tool, called [AmpliconSuite-pipeline](https://github.com/AmpliconSuite/AmpliconSuite-pipeline)**

AmpliconSuite-pipeline wraps the required steps before running AA. Users will enter the AA workflow from different stages. Some will start with fastq files, others will have a BAM file only, and others will
already have the BAM file and CNV seed regions they wish to analyze. We have created this wrapper to allow users to enter 
the workflow from any point. AmpliconSuite-pipeline wraps BWA MEM for alignment, CNVKit for CNV seed identification, and will also invoke
the AA `amplified_intervals.py` filtering script to select/filter/merge CNV seeds which are appropriate for AA.

Ultimately, AmpliconSuite-pipeline can even invoke AA (if installed beforehand), and thus saves users from the hassle associated with preparing everything to run AA on a sample.

If you decide to "got it alone" and not to use AmpliconSuite-pipeline, please carefully read the following points:

#### - Creating the BAM file:
If you are generating your own BAM file, please note that aligners vary in terms of which tags they will add to BAM files. Furthermore,
BAM files hosted on SRA are often stripped of tags needed by AA to correctly identify breakpoints. At this time, we recommend using
**BWA MEM**. However, we have also modified AA so that it will support BAM files created by Illumina's Isaac aligner. Please note that for use with AA, BWA MEM should be run **without** setting the `-m` flag (it's off by default, just don't turn it on).

AmpliconSuite-pipeline will also function on coordinate-sorted CRAM files, [provided that the CRAM reference is in place](https://www.htslib.org/workflow/cram.html).

Please also note that the AA data repo has reference genome fasta files you can align to. We recommend them as they are stripped of alternate contigs, and thus may cause fewer issues with CN estimation in AA.

Don't forget to also create the BAM file index! `samtools index [sample].bam`

#### - Creating the CNV bed file:
It is important to remember that the goal of the CNV calls given to AA are simply to identify locations where focal amplifications may exist - that is, there is a locally increased copy number. This means that CNV callers which sensitively segment the reference genome will perform best for focal amplification detection. CNV callers that undersegment the genome may have good performance in genome-wide tool comparisons, but may not work well with AA.

In the AA publication, ReadDepth was used as the CNV caller for seeding. However, there are much more modern CNV callers available. Some other callers we have used successfully include CNVKit, Canvas, Battenberg, and GATK.

If you generate your own file of CNV calls, *please ensure the estimated copy number of these calls is in the last column of the bed file*. 
Secondly, also ensure that the calls you are using are segmented, and not just raw per-bin estimates. This is not a concern for most users,
however **if you notice there are > 50 CNV seeds going into AA, there may be something wrong**.   

We have found that CNV callers FACETS and ASCAT generally undersegment or create seeds that are too long for AA. **We do not recommend using FACETS or ASCAT for AA seed detection.**

Many modern callers are purity and ploidy aware. Note that while this may help to accurate get globally accurate copy numbers, it can lead to situations where entire chromosome arms surpass your copy number gain threshold for AA! If you are using purity and ploidy corrected calls, please set the `gain` cutoff accordingly.
#

### Filtering and selecting seeds from CNV data
Focal amplifications are somewhat aribtrarily defined as regions of the genome with a CN > 5 and 10 kbp < size < 10 Mbp.

**Selecting appropriate seed regions from the CNV estimates is critical to properly running AA**.

We recommend picking regions which have an estimated CN >= 4.5 and size > 50 kbp, which do not appear amplified due to being parts of repeat elements, and which are not amplified due to karyotypic abnormality. 

AmpliconSuite-pipeline.py calls on multiple filters to ensure that CNV seeds are appropriately selected. Attempting to bypass this filtering and implement some alternative or less rigorous strategy is typically detrimental to getting high-quality focal amplification calls. 

In brief, AmpliconSuite-pipeline identifies seed regions by filtering regions not more than two copies above median chromosome arm ploidy. Filters are also applied based on amplification size and include previously published9 AA filters (amplified_intervals.py) for repetitive and poorly mappable sequences along regions of the reference genome commonly observed to have elevated copy number in non-cancer samples. 

CNV estimates can be imperfect in low-complexity or repetitive regions or appear consistently high when there is karyotypic abnormality. Thus, we have developed modules called `cnv_prefilter.py` and `amplified_intervals.py` to address those issues. They are used by default. If low-complexity/repetive seeds are not filtered from AA, it can cause an exremely long runtime and produce results which are not useful. AA has its own filters for these regions, but it should still be avoided to give them to AA as input.

#

### Constructing unified seed sets from multiple related samples
One common type of analysis is the detection of focal amplifications in samples derived from the same patient or cell line. These may include multiple samples derived from multiregion biopsies, timecourse sequencing on the same individual, or analysis of a cell line before and after drug treatment.

**To maximize the utility of AA on samples from the same source, you should consider merging your AA seed files from related samples.**
We will work on automating this process in the future, but these are the instructions users should follow for now.

#### Option A: `GroupedAnalysisAmpSuite.py` (recommended)
This can be done using the script `GroupedAnalysisAmpSuite.py`. This assumes you are starting from .bam files. Simply place your samples into a three column file, formatted like so 
>`sample_name` `path/to/sample.bam`  `'Tumor' or 'Normal'`

Users can provide two additional columns
> `/path/to/CNV.bed` `/path/to/sample_metadata.json`

These are positional for columns 4 and 5, respectively. If one is ommitted, the `NA` should be placed.

`GroupedAnalysis.py` takes most of the same arguments as `AmpliconSuite-pipeline.py`, but by default `--run_AA` and `--run_AC` will be set.

#### Option B: Manual creation of unified seeds
First run AmpliconSuite-pipeline on each sample **without** setting `--run_AA`. The result for each sample the `[sample]_AA_CNV_SEEDS.bed` file. Next, using `bedtools` or similar, take the relevant seeds files from related samples and merge those bed files. A tutorial for merging bed files can be found [here](http://quinlanlab.org/tutorials/bedtools/bedtools.html) 
(see section entitled "bedtools 'merge'").

Next, place a placeholder copy number estimate into the last column of each row. Since the regions are already filtered, AA does not need to consider the copy number estimate before running its own CN-estimation on each sample.
E.g. `sed -i "s/$/\t999999/" initial_merged_AA_CNV_SEEDS.bed > final_merged_AA_CNV_SEEDS.bed`

The resulting merged file should still end with the suffix `_AA_CNV_SEEDS.bed`, since this suffix has a special meaning in AmpliconSuite-pipeline, and filtering will be skipped (these samples are already filtered.)

You can then run `AmpliconSuite-pipeline.py` with this merged `AA_CNV_SEEDS.bed` file for each of the related samples, now ensuring that each sample is launched on the same collection of regions.

#

### Resource and timing requirements for running AmpliconSuite-pipeline
![image](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/assets/14268531/7f785727-b937-4866-9b42-9ae84c6faee6)

As suggested by the diagram above, which gives timings for typical ~30x coverage samples, you may achieve the best efficiency on your HPC system by breaking up the running of AmpliconSuite-pipeline into stages using the resources necessary for that stage. You can then re-launch the pipeline with different resource allocations using the files generated during the previous steps.

#

### Running AA
Some power-users and other developers may want to run AA on their own, without using the recommended AmpliconSuite-pipeline mode. We assume the user now has a coordinate-sorted BAM file, and a CNV seed BED file (i.e., a BED file of seeds output by AmpliconSuite-pipeline `amplified_intervals.py`).
To check if your BAM file is coordinate-sorted, you can take a peek at the BAM file header by doing 
`samtools view -H your_bamfile.bam | head `. 
Please make sure you know which reference genome it is aligned to so that you can properly specify the `--ref` argument to AA.

For more on running AA, please see the [relevant section of the AA README](https://github.com/virajbdeshpande/AmpliconArchitect#running-ampliconarchitect) or jump down to the [worked example](#worked-example).

Note that AmpliconSuite-pipeline can run AA on its own by setting `--run_AA`, and AA will automatically be called at the end of the preparation process without any additional work by the user. Many AA-specific arguments are exposed in AmpliconSuite-pipeline.py
#

### Interpreting the output
"*How do I know if my sample has ecDNA (or BFB, or other)?*"

**We have recently developed classification methods which take AA output and predict the mechanism(s) of a focal amplification's genesis**.
This method is called **[AmpliconClassifier](https://github.com/jluebeck/AmpliconClassifier)**. 
This tool will output a table describing many important details of the focal amplifications AA discovered. AmpliconClassifier (AC) can be run by setting `--run_AC` when launching AmpliconSuite-pipeline

On its own, AA does not automatically produce a prediction of ecDNA or BFB status. It provides the files though that can be used to make that determination.
For more details about deciphering the AA outputs on your own, please see the [relevant section of the AA README](https://github.com/virajbdeshpande/AmpliconArchitect#outputs).

"*What are the limitations in interpreting the outputs?*"

**Some of the primary limitations to consider with AA go as follows:**

- False negatives: In a panel of 67 FISH-validated amplicons in cancer lines using low-coverage WGS, AA achieved ~87% sensitivity with regards to ecDNA detection (Kim H, et al., Nature Genetics 2020).
- False positives: It is very hard to establish all regions of an arbitrary genome that appear to have bioinformatically circular repeats. We attempt to solve this by providing the similarity score filtering option in AmpliconClassifier, but users should still be wary and thorough with vetting their results.
If they see identical focal amplifications in unrelated samples, this is a red flag.
- Thresholds are unstable: We use by default a threshold of CN=4.5 and size >10kbp when detecting ecDNA. Borderline events near these cutoffs may not be stable in terms of detection.
- ecDNA and BFB are difficult to distinguish bioinformaticall from short read data: While breakage-fusion-bridge cycles are theoretically simple to identify, the presence of additional SVs in, and around the focally amplified region obscures that signature.
Consequently, these events may appear like ecDNA. Furthermore, it has been posited that [BFBs may produce ecDNAs as an intermediate product](https://link.springer.com/article/10.1007/s00412-022-00773-4), meaning that genomically overlapping BFB and ecDNA amplifications may exist, further confounding detection.


#

### Visualizing the output
AA will generate a visualization of the amplicon with edges/CN on its own. The coloring of the edges represents the [orientation of the read pairs](https://github.com/AmpliconSuite/AmpliconArchitect#4-the-sv-view-out_ampliconidpngpdf).

One common task though is, "*I have identified a likely ecDNA in my sample using AA, and would like to visualize it as a cycle*." 

We have developed a tool called **[CycleViz](https://github.com/jluebeck/CycleViz)**, which can produce a circular or linear 
plot of a path or cycle contained in an a AA cycles file. Please note the instructions for running CycleViz in the README, 
as there is a script to reformat the cycles file from AA prior to running CycleViz.
#

### Worked example
Here are some commands for running AA. We'll assume you're starting with a coordinate-sorted BAM file.

- **Run some "preflight checks"**:

```
echo $AA_DATA_REPO
ls $HOME/mosek/
```
None of the above should print an empty string or "no such file" error.

- **Launch AmpliconSuite-pipeline**:

In this specific case, we'll assume you don't have a CNV bed yet, and we'll assume you've installed CNVKit. Please do see the [AmpliconSuite-pipeline README](https://github.com/AmpliconSuite/AmpliconSuite-pipeline) though to check which flags you need to set in your case.
```
AmpliconSuite-pipeline.py -s sample_name  -t number_of_threads --cnvkit_dir /path/to/cnvkit.py --sorted_bam sample.bam --run_AA --run_AC
```

You can still do this even if you already have CNV calls from your own caller, or if you only have the fastq files. Just check the README for the commands to use.

- **Classifying results from multiple runs**:

We can run a classification of multiple AA outputs to determine if ecDNA or other focal amps are present. You will have to install [AmpliconClassifier](https://github.com/AmpliconSuite/AmpliconClassifier). 
```bash
# first make the .input file from multiple runs
make_input.sh /path/to/your/AA_outputs example_collection

# this will create a file called "example_collection.input" 
# next run AC
amplicon_classifier.py --ref GRCh38 --input example_collection.input
``` 

- **Creating a circular visualization of the ecDNA**:

Suppose you find that `amplicon1` of your sample `GBM39` represents an ecDNA+ amplicon, and that in the AA cycles file, the most dominant cyclic path is `Cycle1`. Perhaps you would like to visualize this cycle.

1. Install [CycleViz](https://github.com/jluebeck/CycleViz).
2. Convert your AA cycles file into the correct format 
    
    ```$CV_SRC/convert_cycles_file.py -c GBM39_amplicon1_cycles.txt -g GBM39_amplicon1_graph.txt```
3. Run CycleViz (or LinearViz).
```
# note that we use the converted cycles file
$CV_SRC/CycleViz.py --cycles_file GBM39_amplicon1_BPG_converted_cycles.txt --cycle 1 -g GBM39_amplicon1_graph.txt --ref hg19 --label_segs --gene_subset_file Bushman --gene_fontsize 9
```

At the moment, we do not support adding additional tracks of data into the plot (e.g. RNA-seq, ATAC-seq, etc.), but that is coming soon.
#

### FAQ
- **Can AA detect subclonal ecDNA?**
    - AA is not inherently a subclone-aware method, and will not attempt to compute clone fractions. AA relies on copy number amplifications being strong enough, even among subclonal events, that they can be identified in bulk.
As a result, users should be aware that low copy number amplifications in tumor subclones may be missed by the tool. **Attempting to apply a subclonal CN caller for seeding with AA may result in false positives, long runtimes, and potentially uninterpretable outputs**. Unless every single cell of a tumor is sequenced, there is no way to identify every possible genomic event in that tumor.

    - Users should keep in mind that ultra-rare events below a certain abundance in the tumor may not be currently having a significant effect on the tumor's biology until a selective pressure is applied and they possibly become positively selected for.
If users want to identify those ulta-rare events that may become relevant in the context of future tumor evolution/therapeutic resistance, different tools or technologies should be applied to answer that specific question. 


- **The AA_CNV_SEEDS.bed file was empty, what's wrong?**
    - Likely nothing. If no seed regions are detected, the sample likely has no candidate focal amplifications.
  Do you see any error messages printed to the console or in the log files?


- **Can I use AA with whole-exome sequencing, ATAC-seq, or RNA-sequencing data?**
    - AA is fundamentally incompatible with these data modalities. We only support paired-end WGS data with AA.
    

- **Can I use AA on Circle-Seq data?**
    - AA is not designed for use with targeted sequencing data, and as a result it may crash or produce unreliable amplicon calls. Please proceed very cautiously if deploying AA/interpreting AA outputs on Circle-Seq data.


- **What coverage is needed for AA?**
    - Because the CN of focal amplifications is higher than the background reference, a BAM with 10x coverage will effectively have 50x coverage in a region with CN 10 (assuming 10x coverage for CN=2). Thus, even very low coverage BAM files can be used with AA. Coverage as low as 0.5x still yields robust focal amplification predictions if the amplification is strong and sample purity is high.
  

- **Do I need to use downsampling? Will I detect more SVs with less downsampling?**
  - Within AA, the threshold for number of reads needed to detect a SV scales with effective coverage of the sample. Lower effective coverage (e.g. downsample 10) will have a lower threshold for number of reads needed to identify an SV than higher effective coverage (e.g. downsample 40). For low copy-number focal amplifications (CN < 10), higher effective coverage may perform slightly better due to more sensitive copy number segmentation at high coverage
(consider setting a downsample parameter of 30 or 40). Similarly, ultra-high CN (100+) events may produce cleaner segmentation with default downsampling (10). Additionally, if your sample has a very low properly paired read rate (see AmpliconSuite-pipeline log file), or if it appears that there are many sequencing artifacts in the sample (lawn of brown, teal or pink short vertical lines in the AA image), setting a downsample value of 1 has helped in many instances. While counterintuitive, the lower effective coverage reduces the chance overlapping discordant read pairs from sequencing artifacts get called as a false-positive SV.
    

- **AA has been running for more than 72 hours. What's wrong?**
    - Please ensure that you selected your CNV seeds appropriately. If you use a CN cutoff any lower than 4.5, and size < 10 kbp,
     there can be a great many false positive focal amplifications. Secondly, and very importantly, please ensure you used AmpliconSuite-pipeline to select the CNV seeds appropriately (as opposed to running `AmpliconArchitect.py` on raw, unfiltered copy number data). 
    - If AA is running for > 1 week, please check that you are using the latest version of the tool & update your AA data repo if it is out of date, and consider changing `--cnsize_min 100000 --cngain 5`, and moving up from there with your cutoffs.
     

- **Should copy-number calls provided to AA be corrected for purity and ploidy?**
    - AA's internal analysis of copy number is a coverage-based method that establishes a baseline coverage throughout the genome. 
  Copy number at a given location is also coverage based and reported relative to this baseline coverage. 
  AA itself does not consider purity or ploidy when computing its relative copy number estimates. 
  The goal is not to determine an absolutely correct copy number, but to identify focally amplified regions relative to the rest of the genome.
    - Adding corrections for purity or ploidy to the copy numbers used for seed determination may cause entire chromosome arms to appear above the copy-number cutoff for a focal amplificaiton. We recommend users set a higher CN-cutoff (e.g. `--cngain 8`) if using these kinds of corrected calls.
    - We also suggest that analyzing samples with purity below 30% should be avoided where possible. Only very-strongly amplified focal amplifications will be found. If corrected for purity, any distortions to copy number signal may also become amplified, creating many false-positive seeds.


- **Are there special considerations for FFPE samples?**
    - Recommend trying a few of your FFPE samples with default settings. If you see what looks like a "lawn" of brown edges or a lawn of magenta and teal edges, these are artifactual SVs. In that case, we recommend setting `--AA_insert_sdevs 9.0 --downsample 1`. While somewhat counterintuitive, the benefit of downsampling to 1x is that sequencing artifacts are less likely to stack on top of each other causing AA to predict an artifactual SV.

- **In the amplicon visualizations, what are the different edge colors?**
    - They represent the orientation of the breakpoint. Please see the [relevant section of the AA README](https://github.com/virajbdeshpande/AmpliconArchitect#4-the-sv-view-out_ampliconidpngpdf). Pink and teal are 'inversion-like', brown is 'duplication-like' (jumps backwards in the reference genome), and red is 'deletion-like' (skips forward in the reference genome). Blue corresponds to a 'source edge', which is an SV that has one end mapped in the amplicon and one end at an unknown destination.
  
  
- **What if I use really low CN cutoffs and really small minimum sizes so that I don't miss any ecDNA?**
    - AA will likely stall. A lower CN threshold is also not better, as in that case it may pick up many more regions which represent non-ecDNA events like segmental tandem duplications. When the CN threshold is lowered, the seed sizes also tend to expand greatly - as they may reflect karytoypic events, not focal amplifications. Giving enormous seed regions to AA, (e.g. > 10 Mbp) is strongly not recommended. 


- **How can I construct a data repo for other species?**
  - At this time only hg19, GRCh37, GRCh38, and mm10 are supported. Providing additional support for references of other species requires the construction of an annotation database. 
  Unfortunately that process is quite complicated and requires multiple different annotation files to be available from the UCSC genome browser and other sites. 
  Not all species are as well-annotated as human and mouse so it may not be possible for your species of interest. Construction of a data repo also requires a sizable panel of "normal" WGS samples
  from to use in generating baseline CNV calls throughout the genome. Unfortunately, not all species are as well-annotated as human and mouse so it may not even be feasible. 
  The AA genome annotations are also used for marking low complexity, repetitive regions, low-mappability regions, 
  oncogenes, as well as areas that show high signal across many "normal" samples. Generating such a database takes at minimum a few weeks of dedicated work on our end, without yet factoring in the 
  time it takes to test this new data repo and validate that it works properly. Naturally we must be very judicious about which species we are able to support.
  - The files necessary for the data repo include 
    - A completely resolved reference genome.
    - Coordinates of centromeres.
    - A mappability exclusion file (such as is created by [excluderanges](https://github.com/dozmorovlab/excluderanges)).
    - A database of segmental "super" duplications (large segmental duplications), an annotation created by UCSC.
    - 35bp mappability scores, as generated by the [GEM mappability tool](https://academic.oup.com/bioinformatics/article/36/12/3687/5815974).
    - Genome gene list and coordinates (as .gff)
    - Oncogene list (as .gff, subset of genome gene list)
    - A bed file listing the conserved gain regions found in normal genomes. To generate this, the user needs to take 20 WGS bams from normal samples without ecDNA, call 
    copy number variants against a flat reference, then identify and construct a .bed file of all regions having CN>4.9 in at least 20% of the samples.
  - Following generation of this data repo, then test it on both negative control and positive control samples with validated ecDNA/focal amp status to ensure correctness of the new data repo. 
  - Failure to properly create this database of reference genome annotations can lead to catastrophic levels of false-positives and negatives, rendering the outputs of the tool meaningless.


- **Which viral genomes come with GRCh38_viral?**
  - The oncoviral reference genomes packaged with GRCh38_viral are from the [ViFi publication](https://academic.oup.com/nar/article/46/7/3309/4944397). Broadly these include the papillomavirus sequences in [PaVE](https://pave.niaid.nih.gov/) as well as hepatitis sequences aggregated in [this](https://genome.cshlp.org/content/22/4/593.full) publication. In summary, HPV, HBV, and EBV are present in this reference build.


- **Which HPV16 genome is included with GRCh38_viral?**
    - The HPV16 genome which exists on [PaVE](https://pave.niaid.nih.gov/locus_viewer?seq_id=HPV16REF). For those using NCBI, this is most similar to [AY686584.1](https://www.ncbi.nlm.nih.gov/nuccore/AY686584.1). The PaVE version used in the data repo (which we called hpv16ref_1) differs by 9 point mutations throughout the viral genome, but the sequences are the same length (but do not necessarily share the same exact starting point in that circular genome).


- **How do I run AA on a sample with a viral genome(s) to check for integration/viral genome structure/viral ecDNA?**
  1. If needed, first determine which viral strain(s) are present (e.g. ViFi or [FastViFi](https://github.com/sara-javadzadeh/FastViFi), or other methods).
  2. Align your reads to a human reference fasta file with the relevant viral genomes included.
  3. We provide an oncoviral version of the AA data repo here: https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/. `GRCh38_viral` should the argument given to `--ref`. If using a viral genome not present in the AA data repo, modify the AA data repo fasta file for the corresponding human build to include the viral strain(s) as entries. Re-index fasta file. Recommend copying original data repo to new location to modify, then updating `$AA_DATA_REPO` with the modified version.
  4. The virus will be added to the seed regions if present and amplified. However, if using a viral genome not present in the AA data repo, before launching AA, add the viral genome to your `AA_CNV_SEEDS.bed` file.


- **How do I tell if an amplification is due to ecDNA or segmental tandem duplication?**
    - For low CN < 4.5, it's difficult to tell with short reads. However, keep the following in mind - when segmental tandem duplications occur, the same breakpoints must be reused every time it is repeated. When a "cyclic" structure accumulates
    in high CN, this would involve the exact reuse of the same breakpoints many, many times in the case of segmental tandem dups. The simpler hypothesis as CN increases, is that it is mediated through an extrachromosomal DNA mechanism. 
   
 
- **Can AA determine the difference between HSR and ecDNA? Can it find integration points?**
    - AA itself does not determine if an ecDNA amplification exists in double-minute form or HSR (homogeneously staining region) form. From previous observations, ecDNA maintains much of its structure when it integrates into the genome (Turner 2017, *Nature*, Deshpande 2019, *Nat. Comms.*, Luebeck 2020 *Nat. Comms.*). While it must break the circle to linearize, the SVs that originally formed the circle are still there (the original circle-forming SV is not 're-used'). Unfortunately without some sort of imaging data (FISH) or long-range sequencing (Bionano, PacBio, Nanopore), it is not possible to reliably make that determination from AA. In cancer cell lines, such as COLO320DM/COLO320HSR there may be no ecDNA in the HSR version of the cell line, however in real tumors which contain heterogeneous populations of HSR-integrated and double minute (DM) forms of ecDNA, there is likely a reservoir of cells carrying the DM-form even if HSR appears to be dominant, and under the appropriate conditions the DM form may become dominant again (Nathanson 2014, *Science*) provides a nice illustration of this). With this in mind, determining if ecDNA is in DM form or also is HSR-integrated is slightly less materially important, however it will be important to one day have computational methods capable of tackling this task.


- **There's a region I want AA to examine, but it didn't appear in the CNV seeds, what do I do?**
    - You can force AA to run on a region by creating a BED file with that region specifically included, and then invoking AA directly (do not use AmpliconSuite-pipeline or `amplified_intervals.py`). This can also be used to standardize the same AA run on multiple different samples.
 

- **In the cycles file, what do the paths that begin/end with '0+' mean?** 
    - These indicate that the path is non-cylic, and proceeds or is preceeded by the next reference coordinate.

#
 
