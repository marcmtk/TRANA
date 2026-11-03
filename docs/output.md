# genomic-medicine-sweden/TRANA: Output

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
