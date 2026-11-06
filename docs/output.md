# marcmtk/TRANA: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Barbell](#barbell) - Demultiplexing and trimming
- [Deacon](#deacon) - Human read removal
- [EMU call check](#emu-call-check) - Calls with and without reads that match their references almost exactly
- [Savont and classifier comparison](#savont-and-classifier-comparison) - ASVs and their classification, side by side with EMU
- [Read buckets](#read-buckets) - Where every read goes, and what the reads behind each EMU call are made of
- [FastQC](#fastqc) - Raw read QC
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Barbell

<details markdown="1">
<summary>Output files</summary>

- `barbell/`
  - `barbell_read_counts.tsv`: reads per barcode, with the sample it belongs to and a final `run_input` row with the reads in the run. Barcodes without a sample act as blank wells.
  - `samples_without_reads.txt`: samples in the sheet that received no reads. These aren't analysed, and the pipeline also warns about them.
  - `barbell_log.log`: barbell's log, including the patterns found.
  - `barbell_{annotation,filtered,pattern_per_read}.tsv.gz`: per-read tables, with `--save_barbell_annotation`.
- `fastq_pass_merged/fastq_pass_demuxed/`: reads per sample, with `--save_merged_reads`. These still contain any human reads.

</details>

### Deacon

<details markdown="1">
<summary>Output files</summary>

- `deacon/`
  - `host_removal_summary.tsv`: reads in and reads removed per sample.
  - `*_dehosted.json`: deacon's summary per sample.
  - `*_dehosted.fq.gz`: the remaining reads, with `--save_dehosted_reads`.

</details>

Samples with no reads left aren't analysed, and the pipeline warns about them.

### EMU call check

<details markdown="1">
<summary>Output files</summary>

- `results/emu_call_check.tsv`: with `--check_emu_calls`, one row per sample and EMU call.
  - `reads_assigned`: reads EMU assigns to the species with probability at least `--check_min_probability` (0.5).
  - `near_exact_reads`: those of them whose best alignment to a reference of the species scores at least 1 - `--check_max_error` (0.99) times the read length, where score = matches - mismatches - inserted bases - deleted bases.
  - `status`: `confirmed` with at least `--check_min_reads` (2) near-exact reads, otherwise `unconfirmed`. `not_assessed` for every call in a sample where less than `--check_min_share` (0.2) of the main species' reads are near-exact, as read quality is then too low for the check.
  - `calibration`: `uncalibrated` for now. The thresholds were chosen on one run (2026-07-02) and haven't been checked on another run with known composition, so use the labels as a pointer, not as a result. The master report marks the card the same way.
  - `main_species_near_exact_share`: that share, for the call with the most assigned reads.
- The master report lists each sample's unconfirmed and confirmed calls, and `results/classifier_comparison.tsv` has the label in `emu_call_check`, which is empty without `--check_emu_calls`.
- EMU's alignments stay in its work folder as `*_emu_alignments.bam`. They are only published, as SAM, with `--keep_files`.

</details>

An unconfirmed call has no read that fits one of its references well. On the run used to evaluate the check, these were calls made from chimeric reads, but a species with no reference within about 1% would look the same. A confirmed call can still be a contaminant or a species that 16S can't separate from a relative, such as _E. fergusonii_ next to _E. coli_. See `docs/evaluation/emu_chimera_step.md`.

### Savont and classifier comparison

<details markdown="1">
<summary>Output files</summary>

- `results/classifier_comparison.tsv`: one row per sample and species, with the abundance from each classifier.
  - `detected_by`: `both`, `emu` or `savont`.
  - `emu_call_check`: the label from the [EMU call check](#emu-call-check), with `--check_emu_calls`.
  - `savont_tied_genera`: set when references from several genera tie at the lowest edit distance for an ASV behind that species. A pair that 16S can't separate, such as _E. coli_/_Shigella_, will show up here, but so will a mislabelled reference, so treat any other pair as unconfirmed at genus level.
- `results/savont/`
  - `*.species_abundance.tsv`, `*.genus_abundance.tsv`: relative abundance.
  - `*.asv_mappings.tsv`: every reference tied at the lowest edit distance for each ASV.
  - `*.final_asvs.fasta`, `*.feature-table.tsv`: ASV sequences and their read depth.
  - `*_savont_asv.log`, `*_savont_classify.log`: logs.
- `pipeline_info/emu_db_summary.tsv`: number of references and checksums of the database that was used.

</details>

### Read buckets

<details markdown="1">
<summary>Output files</summary>

- `read_buckets/`
  - `read_buckets.html`: the report, linked from the master report.
  - `read_buckets_summary.tsv`: one row per sample. Reads by fate (`deacon`, `short`, `long`, `qfail`, `filter_other`, `emu_not_sampled`, `emu_not_assigned`, `emu_assigned`), reads by bucket for all reads (`bucket_*`) and for the reads EMU got (`emu_bucket_*`), and reads by starting strand.
  - `read_buckets_calls.tsv`: one row per sample and EMU call. EMU's assignment probabilities summed by bucket (`sum_p_*`) and reads assigned with probability at least 0.5 by bucket (`reads_p50_*`), the call check label with `--check_emu_calls`, and `sample_sheet` (`expected` or `not expected`) when the sheet has an `expected` column.
- `read_buckets/samples/`
  - `*.read_buckets.tsv.gz`: one row per read barbell gave the sample: length, mean accuracy, what host removal, filtlong and EMU did with it, its structure and bucket.
  - `*.read_buckets.json`: the references used, each EMU call by bucket, and whether filtlong's choices matched the length and accuracy rule.
- MultiQC gets two bar charts: the fate of each sample's reads, and the buckets of the reads EMU got.

</details>

The report follows every read from `fastq_pass` through barbell, deacon, filtlong and EMU, and checks that the counts add up at every stage. A stage that doesn't add up, or a read filtlong treats differently from its own length and accuracy rule, is listed at the top of the report.

With barbell and `--barbell_kit SQK-MAB114-24-16S`, it also aligns every read to the sample's own 16S sequences (savont's ASVs, named after their closest database references) and to database references of each species EMU called in the sample (one if an ASV has the species, up to three if not), and names the read's structure: full length, truncated, a fold-back onto the other strand, a same-strand repeat, an unaligned tail, or a chimera, where two species with one breakpoint explain the read clearly better than any one species does. Positions are given on _E. coli_ 16S. Barcode patterns are written in the notation of Beeloo et al. (barbell), extended for the 16S insert.

The chimera call was calibrated on one run (2026-07-02), with pure cultures as negatives and reads of two cultures spliced together as positives, and with the consensuses of that run's pure cultures as references for every sample. With the references above instead, 96.8–99.9% of each sample's reads in that run landed in the same bucket, no pure-culture read was called chimeric (against 2–6 per culture), and the same nine unexpected EMU calls rested at least half on chimeric reads. Chimeras between close relatives are under-counted: in a mix with _K. pneumoniae_ and _K. aerogenes_, 379 chimeric reads were found against 1,198 with the culture references, and the other four mixes had between 11% fewer and 0.3% more. A pure culture is only a weak negative control here, because its references are those of its own species. Treat the rates as lower bounds. The chimera call is uncalibrated until it has been checked on another run with known composition, and the report says so.

On the 2026-07-02 run (38,000–59,000 reads per sample), the structure step took 0.6–4 CPU minutes per pure culture and 2–29 per mix, 78 CPU minutes for the whole run, against 7.6 CPU hours of EMU for one of the mixes alone. Time grows with the number of references, so with mixes that have many species. `--skip_read_buckets` turns the report off.

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
