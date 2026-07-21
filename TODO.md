# AmpliconSuite-pipeline TODO

## Make standalone-installed commands discoverable

The README shows bare commands such as `AmpliconSuite-pipeline.py`, but a full
`install.sh` run currently installs dependencies, clones AA and AC, and exports
`AA_SRC`/`AC_SRC` without placing PAA, GroupedAnalysis, AA, or AC on `PATH`.
Internal PAA calls still work because they use the source variables, but direct
commands fail after users leave the repository directory.

- Choose a user-level installation mechanism: an idempotent managed `PATH`
  block or launchers in a user bin directory. Do not use the current
  `setup.py` as-is; its distribution name/version metadata is incomplete.
- Wire `AmpliconSuite-pipeline.py`, `GroupedAnalysisAmpSuite.py`,
  `AmpliconArchitect.py`, and `amplicon_classifier.py`; verify the pip-installed
  `BFBArchitect.py` entry point too.
- Support both Bash and Zsh, preserve existing `AA_SRC`/`AC_SRC` selections,
  quote paths, and make `--uninstall` remove only installer-managed entries.
- Keep `--finalize_only` safe for Conda and container users, where the downloaded
  installer may not reside in the PAA repository.
- Add an isolated-`HOME` test that opens a fresh shell and checks `command -v`
  plus `--version` for each promised command. Reconcile README examples with
  the behavior actually tested.

## Add fast BAM QC and automatic AA tuning

Replace the default whole-file `samtools flagstat` pass with a faster,
reproducible QC workflow while retaining a full-scan fallback.

- Run a fast BAM/CRAM integrity check (for example, `samtools quickcheck`) for
  headers and truncation. Keep this separate from statistical sampling;
  sampling cannot prove that the entire file is readable.
- With the coordinate index already created by PAA, sample deterministic random
  genomic windows, similar to AA's coverage-statistics sampling. Collect enough
  primary reads to estimate the properly-paired rate and report sample size and
  uncertainty. Define the denominator explicitly and validate estimates against
  `samtools flagstat` on representative BAMs.
- Measure foldback-like pair orientations and short discordant-pair artifacts.
  Distinguish widespread library artifacts from localized biological foldback
  signal before changing parameters.
- Translate validated QC findings into conservative AA settings. Preserve
  explicit user arguments; candidate automatic responses include raising
  `--AA_insert_sdevs`, raising `--foldback_pair_support_min`, and, for strongly
  artifacted libraries, reducing AA downsampling coverage.
- Log sampled metrics, thresholds, resulting parameter changes, and their
  provenance in a machine-readable QC output as well as the PAA log.
- Test indexed BAM and CRAM inputs, truncated files, sparse/targeted data,
  low-proper-pair libraries, artifacted libraries, and normal WGS. Exercise both
  standalone and GroupedAnalysis paths and benchmark runtime against flagstat.
