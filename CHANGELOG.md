# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Version 2.0.0 of the marcmtk hard fork, based on genomic-medicine-sweden/TRANA v1.1.0.

### Added

- Demultiplex and trim `--merge_fastq_pass` runs with the marcmtk fork of barbell (`--skip_barbell` to opt out)
- Remove human reads with deacon, using the nf-core module (`--deacon_index`, `--skip_host_removal`)
- Run the marcmtk fork of savont on the same reads as EMU (`--run_savont`)
- `results/classifier_comparison.tsv` with both classifiers' calls and cross-genus ties in savont
- `--check_emu_calls`, off by default, which labels EMU calls without near-exact reads as unconfirmed in `results/emu_call_check.tsv`, the master report and the classifier comparison. EMU's alignments are then kept as BAM
- `PREPARE_EMU_DB`, which accepts a folder or `.tar.gz`, pads `taxonomy.tsv` for savont and drops the references in `--db_exclude_ids`
- `workstation` profile and `pixi.toml` for running with podman on one machine
- Warnings when a sample is dropped for having no reads
- `SAVONT_INDEX`, which builds savont's index once from the prepared database, so that `SAVONT_CLASSIFY` no longer writes beside the shared database. That write could crash concurrent tasks and made `-resume` rerun emu
- The master report links run-level results, savont, deacon, barbell, filtlong and the processed NanoPlot reports, showing only what the run produced
- A read bucket report (`read_buckets/read_buckets.html`) that follows every read through barbell, deacon, filtlong and EMU, names each read's structure on 16S (full length, truncated, fold-back, repeat, chimera) with barbell's MAB 16S kit, and shows what the reads behind each EMU call are made of. Optional `sample_type` and `expected` columns in the barcodes sample sheet group the samples and mark calls as expected or not. `--skip_read_buckets` turns it off

### Changed

- The bundled database is now the 2026 EMU database (rrnDB v5.10 and NCBI RefSeq 16S, built by @UFDuttonLab)
- `--db_exclude_ids` is on by default with the audit list for that database (`--skip_db_exclusion` to opt out)
- Krona plots use `ktImportText` with the lineage from the EMU report. `--krona_taxonomy_tab` and the bundled Krona taxonomy are gone, as 2,821 of the database's 27,670 taxids were missing from it
- EMU 3.5.4 to 3.6.2
- The database is staged as a process input rather than read from a host path
- `--save_merged_reads` is off by default, as those reads still contain any human reads
- Nextflow `>=25.04.0`, as the deacon module uses topic channels
- The pipeline is now `marcmtk/TRANA` 2.0.0, a hard fork maintained separately from genomic-medicine-sweden/TRANA. CI runs on the fork, without conda since barbell and savont aren't on bioconda, and the tool images are pinned by digest

## [v1.1.0]

### Changed

- Made the `lr:hq` minimap2 profile a valid option to the `--seqtype` parameter

## [v1.0.0]

### Added

- Added `generate_yaml` process to trana workflow
- Added outdir to workflows for `generate_yaml`

### Fixed

- Replaced `check_max()` with `process.resourceLimits`
- Added stub blocks to `TRANSLATE_TAXIDS`, `ASSIGNMENT_HEATMAP`, and `EMU_COMBINE_OUTPUTS` for CI stub runs

### Changed

- Downsampled reads are now included in the second round of Nanoplot
- Fixed broken container image URL for NanoPlot
- Make parameters_schema configurable for nf wrappers
- Replace `sample_size` with `downsample_n_reads` to avoid ambiguity

## [v0.6.0]

### Added

- Added diff between local and remote modules

### Fixed

- Fixed typo in readme
- Fixed bioconda channel for generate_master
- Changed from "TACO" to "TRANA" in the metro map in the README.

### Changed

- Changed `seqtk` out put publishing to be optional via `save_downsampled_reads` param
- Updated nanoplot to 1.46.1 module from nf-core and modified this to ouput with prefix
- Updated cond env files
- Updated contributors in nextflow.config and in ro-crate-metadata.json
- Updated Emu to 3.5.4 to fix [this float diffing issue](https://github.com/treangenlab/emu/issues/90) causing some datasets to crash.

## [v0.5.0]

### Added

### Fixed

- Fixed `filtlong` bug that crashes pipeline when all reads are filtered out with `.branch` operator
- Fixed fastqc data missing in multiqc report
- No v in version in nextflow.config

### Changed

## [v0.4.1]

### Added

### Fixed

- Fixed version in nextflow.config

### Changed

## [v0.4.0]

### Added

## [v0.4.0]

### Added

- Added a `generate_master_html` python script that creates `master.html` file
  containing a table of samples with corresponding pointers to each html output
  file
- Added repective `GENERATE_MASTER_HTML` process
- Added `cmd.config`
- Added `params.trace_timestamp` to `nextflow.config`
- Added `changelog_update_reminder` GA workflow
- Added optional ability to save merged reads
- Added capability to emit nanoplot nanostats text file
- Added clear commenting for modules and subworkflows
- Added stubs for all local modules
- Added `when` operator for nanoplot process
- Added workflow that publishes docker image for python-related processes
- Added contributors
- Added versioning to `merge_barcodes_samplesheet.py`
- Added output argument to `generate_master_html.py`
- Added prefix to `master.html` filename
- Added `seqtk_sample` process
- Added `Makefile` rules for `precommit`, `lint`, `schema`, `test`,
  `test-cli-fastq`, `test-cli-samplesheet`, as well as a `check` command that
  runs all checks
- Added `nf-test`
- Added `pre-commit` checks, which includes the `prettier` tool
- Added GitHub Actions workflows for `nf-core` linting checks
- Updated nf-core modules to include nf-tests and added minimal test
- Added emu-combine-output
- Removed authors-field in nextflow.config.
- Added a module to output the result as a phyloseq object
- Added figure for pipeline
- Added an nf-test for short-read, paired-end data
- Added tests for the combinations of adapter trimming and quality filtering parameters.
- Added test for downsampling
- Added running of `nextflow lint` to CI
- Added module for heatmap for emu likelihood file
- Added barplots for control comparison
- Added module for translation of tax-ids into scientific name for heatmap generation
- Logs for local modules
- result directory for seqtk
- Added a config that can be used for SGE

### Fixed

- Readded `NANOPLOT2` as `NANOPLOT_PROCESSED_READS`
- Fixed Dockerfile context
- Moved nanostats_unprocessed process execution into seqtype SR if statement
- Conditionally emit nanostats unprocessed/processed to avoid undefined output error when using --seqtype SR
- Fixed a broken configuration file in `modules/local/emu/abundance/meta.yml`
- Fixed failing pulls of some singularity containers by settings singularity.cacheDir
- Fixed broken nanoplot code. Readded prefixes to the module.
- Fixed broken `--help` screen by adding missing nf-schema plugin and validation config.
- Fixed broken e-mail sending, where a (corrupt) e-mail would be sent at the start of a pipeline rather than at the end.
- Fixed a pipeline ordering error in the main workflow where Nanoplot outputs were accessed before executing the module.
- Fixed various broken process calls in the main Trana workflow found when adding tests for parameter combinations.
- Updated schema
- Fixed input channel for seqtk
- Fixed names for containers
- Fixed seqtk links

### Changed

- `merge_barcodes_samplesheet.py` can now handle custom barcodes.
- Provided option to `save_merged_reads`
- Provided ability to overwrite files in `publishDir` (`params.publish_dir_overwrite`)
- Cleaned up input preprocessing steps
- Updated `master_template.html` for nanoplot prefix and output changes
- Changed `merge_barcodes_samplesheet.py` container to nf-core (temporary)
- Emit `master.html`
- Changed project name to TRANA (`genomic-medicine-sweden/trana`)
- Updated nf-core template to version 3.2.0
- Updated the `ci.yml` config to run test via `nf-test` in addition to the normal test run via cli
- Updated `README.md` with info about running `make install` to gunzip gzipped
  assets and useful `make` commands for developers.
- Updated workflow image
- Made phyloseq an option to protect against errors when using other databases than emu´s premade.
- Added peoples github names to README for credit
- Added nanoplot unprocessed reads and processed reads to multiqc. Short info about this in README.
- Added steps in `PIPELINE_INITIALISATION` from nf-core template, and replaced workflow.onComplete code
  with call to `PIPELINE_COMPLETION` subworkflow
- Updated images for github
- Added versions for local scripts
- Changed `PIPELINE_INITIALISATION` subworkflow to be modularised
- Changed `bin/merge_barcodes_samplesheet.py` to take csv as input instead of tsv as stated in `README` (safer option)
- Name of pipeline was change to TRANA

## [v0.1.0]

### Added

- Support for Nanopore reads. Reads can be in barcode directories or premerged into samplespecific fastq-files.
- Krona
- emu
- FiltLong
- Porchop_abi
- NanoPlot
- Multiqc
- fastqc
- Adaptation to nextflow
- All tools are containerised with singularity
