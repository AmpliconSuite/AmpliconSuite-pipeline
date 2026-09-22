# CRAM inputs and reference caches

Use `--bam sample.cram` for standalone analysis, or supply CRAM paths in the
alignment-file column of the grouped input file. CRAM 3.1 and higher require
samtools >=1.21 in AmpliconSuite. Preflight checks the file's format version and
external samtools version before downloading references or decoding an alignment.
Future CRAM formats may require additional decoder updates; this minimum is not
a guarantee of future-format support. CNVKit uses its own pysam/HTSlib decoder;
changing external samtools does not update pysam.

Both entry points use the shared CRAM preflight. Grouped analysis checks distinct
tumor and normal CRAM inputs before starting sample jobs. Reference lookup failures
stop the run and print a reference-fetch command. The fast decoding check reads
one alignment without requiring an index; it cannot guarantee that every reference
needed later is available. Failed indexing also stops the pipeline.

## Populate a persistent cache

Run this with the AmpliconSuite Python environment, from the repository directory
(or from any directory if `paalib` is installed in that environment):

```bash
export REF_CACHE="$HOME/.cache/hts-ref/%2s/%2s/%s"
export REF_PATH="$REF_CACHE"
python -m paalib.cram_compat --fetch_cram_refs /path/to/sample.cram
```

You can give multiple CRAM paths. This command uses the same downloader as the
standalone pipeline's `--fetch_cram_refs`: reference MD5s come from the CRAM header,
missing sequences are downloaded from ENA, and downloaded bytes are checksum
validated before installation. Existing files are skipped. Subsequent runs reuse
the cache; there is no reason to download it afresh for every analysis. Keep the
same environment settings when starting standalone or grouped analysis. Omit
`--fetch_cram_refs` to perform only the preflight.

The reference must match the sequence used when encoding the CRAM, including any
needed alternative/decoy sequences. An AA data repository FASTA alone need not
contain all of them. Paths recorded in CRAM headers may point to a different
machine; a checksum-addressed cache avoids dependence on those paths.

## Docker and Singularity launchers

`docker/run_paa_docker.py`, `singularity/run_paa_singularity.py`, and
`singularity/run_ga_singularity.py` automatically prepare references when any
tumor or normal input is a CRAM. They use the host `REF_CACHE` template if set,
otherwise `$HOME/.cache/hts-ref/%2s/%2s/%s`. Additional local checksum templates
from `REF_PATH` are mounted read-only and reused. Remote `REF_PATH` URLs and FASTA
files are not scanned by this cache-preparation step; missing checksum entries
are downloaded from ENA into `REF_CACHE`.

Preparation runs inside the selected image using the shared reference fetcher;
the host does not need samtools or pysam. It checks the CRAM header's reference
checksums and downloads only absent entries. It then checks decoding and mounts
the reference stores read-only for analysis, with container-side `REF_CACHE` and
`REF_PATH` set explicitly. Preparation failures stop the launcher before analysis.
Existing cache files are not checksum-scanned or repaired.

Docker CRAM runs automatically use the host UID/GID so both preparation and
analysis can access the same private cache. Singularity retains the launching
user's identity. The cache persists independently of output cleanup. No separate
fetch flag or manual cache mount is needed with these launchers:

```bash
export REF_CACHE="$HOME/.cache/hts-ref/%2s/%2s/%s"  # optional: this is the default
python docker/run_paa_docker.py --image YOUR_IMAGE \
  -s sample -t 4 --bam sample.cram --ref GRCh38 -o output
python singularity/run_paa_singularity.py --sif YOUR_IMAGE.sif \
  -s sample -t 4 --bam sample.cram --ref GRCh38 -o output
python singularity/run_ga_singularity.py --sif YOUR_IMAGE.sif \
  -i samples.txt -t 4 --ref GRCh38 -o grouped_output
```

Docker's optional `--image` selects the same image for preparation and analysis;
the default remains `jluebeck/ampliconsuite-pipeline`. Use updated images containing
this reference-fetch command and samtools >=1.21. Build-time checks enforce the
minimum for new Docker and Singularity images; older images must be rebuilt.

## Direct container commands

Mount the cache explicitly and set **container-side** `REF_CACHE` and `REF_PATH`.
Do not assume the host home directory or host environment is available. In
particular, the provided Singularity launchers use `--no-home --cleanenv` and now
configure the cache mounts explicitly. The direct commands below show how to provide
these mounts yourself without a launcher. Use the same mounts/environment when
invoking either pipeline entry point, alongside its usual input/output/data mounts.

These examples require an image containing the updated CRAM module. Set
`IMAGE` or `SIF` to your rebuilt image. Check `samtools --version` inside it:
current build definitions specify 1.21, but older images can still contain 1.13.
The repository Conda recipe requires `samtools >=1.21` and `pysam >=0.23.3`.

Docker, downloading to a persistent host directory with the host user's ownership:

```bash
CACHE_DIR="$HOME/.cache/hts-ref"
mkdir -p "$CACHE_DIR"
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$CACHE_DIR:/cache" -v /absolute/input/directory:/input:ro \
  -e 'REF_CACHE=/cache/%2s/%2s/%s' -e 'REF_PATH=/cache/%2s/%2s/%s' \
  -w /home/programs/AmpliconSuite-pipeline-master \
  --entrypoint python3 "$IMAGE" \
  -m paalib.cram_compat --fetch_cram_refs /input/sample.cram
```

Singularity/Apptainer, also writing directly to the persistent host directory:

```bash
CACHE_DIR="$HOME/.cache/hts-ref"
mkdir -p "$CACHE_DIR"
singularity exec --cleanenv --no-home \
  --bind "$CACHE_DIR:/cache" --bind /absolute/input/directory:/input:ro \
  --env 'REF_CACHE=/cache/%2s/%2s/%s' --env 'REF_PATH=/cache/%2s/%2s/%s' \
  --pwd /home/programs/AmpliconSuite-pipeline-master "$SIF" \
  python3 -m paalib.cram_compat --fetch_cram_refs /input/sample.cram
```

Once populated, mount `/cache` read-only for analysis (`:/cache:ro`). Ensure the
container user can read the files; for downloads it must also be able to write
the host cache. Docker's `--user` setting above avoids creating cache files owned
by another user. Keep this cache separate from disposable container filesystems.
