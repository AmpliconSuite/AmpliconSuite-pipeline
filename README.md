## PrepareAA

A quickstart tool for [AmpliconArchitect](https://github.com/virajbdeshpande/AmpliconArchitect). Performs all preliminary steps (alignment, CNV calling, seed interval detection) required prior to running AmpliconArchitect.

Requires the following tools to be installed beforehand:
- bwa mem
- samtools
- [freebayes](https://github.com/ekg/freebayes)
- [Canvas](https://github.com/Illumina/canvas)
- AmpliconArchitect

Files in hg19/ folder must be placed into $AA_DATA_REPO/hg19 folder prior to using.

Canvas can be tricky to install. Unfortunately PrepareAA only supports hg19 currently, as the standalone version of Canvas currently only seems to support hg19.
