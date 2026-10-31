/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                       } from 'plugin/nf-schema'
include { softwareVersionsToYAML                 } from '../subworkflows/nf-core/utils_nfcore_pipeline/main.nf'
include { paramsSummaryMultiqc                   } from '../subworkflows/nf-core/utils_nfcore_pipeline/main.nf'
include { methodsDescriptionText                 } from '../subworkflows/local/utils_nfcore_trana_pipeline/main.nf'
include { GENERATE_YAML                          } from '../modules/local/generate_yaml/main.nf'
include { GENERATE_MASTER_HTML                   } from '../modules/local/generate_master_html/main.nf'
include { EMU_ABUNDANCE                          } from '../modules/local/emu/abundance/main.nf'
include { KRONA_KTIMPORTTEXT                     } from '../modules/nf-core/krona/ktimporttext/main.nf'
include { EMU_REPORT_TO_KRONA                    } from '../modules/local/emu/report_to_krona/main.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS            } from '../modules/nf-core/custom/dumpsoftwareversions/main.nf'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main.nf'
include { FASTQC                                 } from '../modules/nf-core/fastqc/main.nf'
include { CUTADAPT                               } from '../modules/nf-core/cutadapt/main.nf'
include { NANOPLOT as NANOPLOT_UNPROCESSED_READS } from '../modules/nf-core/nanoplot/main.nf'
include { NANOPLOT as NANOPLOT_PROCESSED_READS   } from '../modules/nf-core/nanoplot/main.nf'
include { PORECHOP_ABI                           } from '../modules/nf-core/porechop/abi/main.nf'
include { SEQTK_SAMPLE                           } from '../modules/nf-core/seqtk/sample/main.nf'
include { FILTLONG                               } from '../modules/nf-core/filtlong/main.nf'
include { EMU_COMBINE_OUTPUTS                    } from '../modules/local/emu/combine_outputs/main.nf'
include { COMBINE_REPORTS                        } from '../modules/local/phyloseq/main.nf'
include { PHYLOSEQ_OBJECT                        } from '../modules/local/phyloseq/main.nf'
include { ASSIGNMENT_HEATMAP                     } from '../modules/local/assignment_heatmap/main.nf'
include { CTRL_COMPARISON                        } from '../modules/local/ctrl_comparison/main.nf'
include { TRANSLATE_TAXIDS                       } from '../modules/local/translate_taxids/main.nf'
include { PREPARE_EMU_DB                         } from '../modules/local/prepare_emu_db/main.nf'
include { DEACON_FILTER                          } from '../modules/nf-core/deacon/filter/main.nf'
include { SAVONT_ASV                             } from '../modules/local/savont/asv/main.nf'
include { SAVONT_CLASSIFY                        } from '../modules/local/savont/classify/main.nf'
include { COMPARE_CLASSIFIERS                    } from '../modules/local/compare_classifiers/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TRANA {

    take:
    ch_samplesheet  // channel: samplesheet read in from --input
    ch_reads        // channel: reads from PIPELINE_INITIALISATION
    outdir

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    // MODULE: Unpack the emu database into the form both emu and savont read
    PREPARE_EMU_DB(
        file(params.db, checkIfExists: true),
        params.db_exclude_ids && !params.skip_db_exclusion ? file(params.db_exclude_ids, checkIfExists: true) : []
    )
    ch_versions = ch_versions.mix(PREPARE_EMU_DB.out.versions)
    ch_emu_db = PREPARE_EMU_DB.out.db

    // MODULE: Run FastQC
    FASTQC (ch_reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // MODULE: Remove human reads with deacon
    if (params.skip_host_removal) {
        ch_reads.set { ch_dehosted_reads }
    } else {
        if (!params.deacon_index) {
            error "Human read removal needs --deacon_index, e.g. the panhuman-1 index from 'deacon index fetch panhuman-1'. Set --skip_host_removal to run without it."
        }
        DEACON_FILTER(
            ch_reads.map { meta, reads -> [ meta, file(params.deacon_index, checkIfExists: true), reads ] }
        )
        ch_versions = ch_versions.mix(
            DEACON_FILTER.out.versions_deacon
                .first()
                .map { process, tool, version -> "\"${process}\":\n    ${tool}: ${version}\n" }
                .collectFile(name: 'deacon_versions.yml')
        )
        DEACON_FILTER.out.log
            .map { meta, summary ->
                def stats = new groovy.json.JsonSlurper().parse(summary)
                [ meta.id, stats.seqs_in, stats.seqs_removed, stats.seqs_removed_proportion ].join('\t')
            }
            .collectFile(
                name: 'host_removal_summary.tsv',
                seed: 'sample\treads_in\treads_removed\tproportion_removed',
                newLine: true,
                sort: true,
                storeDir: "${params.outdir}/deacon"
            )

        // Samples with nothing left can't be classified, and must not vanish unnoticed
        DEACON_FILTER.out.fastq_filtered
            .branch { meta, reads ->
                passed: hasReads(reads)
                failed: true
            }
            .set { ch_deacon_filtered }
        ch_deacon_filtered.failed
            .subscribe { meta, reads -> log.warn "No reads left in sample ${meta.id} after human read removal, so it won't be analysed" }
        ch_deacon_filtered.passed.set { ch_dehosted_reads }
    }

    if (params.seqtype in ["map-ont", "lr:hq"]) {

        // MODULE: Run NANOPLOT_UNPROCESSED_READS
        NANOPLOT_UNPROCESSED_READS(ch_reads)
        ch_versions = ch_versions.mix(NANOPLOT_UNPROCESSED_READS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            NANOPLOT_UNPROCESSED_READS.out.txt.collect{ it[1] }
        )

        // Optional adapter trimming
        if (params.adapter_trimming) {
            // MODULE: Run Porechop to trim ONT adapters
            PORECHOP_ABI(ch_dehosted_reads, [])

            PORECHOP_ABI.out.reads.map {
                meta, reads -> [meta + [single_end: 1], reads]
            }.set{ ch_optionally_trimmed_reads }

            ch_versions = ch_versions.mix(PORECHOP_ABI.out.versions.first())
            ch_multiqc_files = ch_multiqc_files.mix(
                PORECHOP_ABI.out.log.collect{ it[1] }
            )
        } else {
            ch_dehosted_reads.set{ ch_optionally_trimmed_reads }
        }

        // Optional quality filtering
        if (params.quality_filtering) {
            // MODULE: Run filtlong to filter on read length
            FILTLONG(ch_optionally_trimmed_reads.map {
                meta, reads -> [meta, [], reads]
            }).reads.set { ch_processed_reads }

            // Split channels into read files with or without read content
            ch_processed_reads
            .branch { meta, reads ->
                passed: hasReads(reads)
                failed: true
            }.set { ch_processed_filtered_reads }

            ch_processed_filtered_reads.failed
                .subscribe { meta, reads -> log.warn "No reads left in sample ${meta.id} after quality filtering, so it won't be analysed" }
            ch_processed_filtered_reads.passed.set { ch_processed_passed_reads }

            ch_versions = ch_versions.mix(FILTLONG.out.versions.first())
            ch_multiqc_files = ch_multiqc_files.mix(
                FILTLONG.out.log.collect{ it[1] }
            )
        } else {
            ch_optionally_trimmed_reads.set{ ch_processed_passed_reads }
        }

    } else if (params.seqtype == "sr") {

        // MODULE: Run cutadapt to trim short-reads adapters
        if (!params.skip_cutadapt) {
            CUTADAPT(ch_dehosted_reads)
            CUTADAPT.out.reads.set { ch_processed_passed_reads }
            ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())
            ch_multiqc_files = ch_multiqc_files.mix(
                CUTADAPT.out.log.collect{ it[1] }
            )
        } else {
            ch_dehosted_reads.set { ch_processed_passed_reads }
        }
    } else {
        error "Invalid seqtype. Please specify either 'map-ont', 'lr:hq' or 'sr'."
    }

    if ( params.downsample_n_reads ) {
        // MODULE: Downsample reads
        SEQTK_SAMPLE(
            ch_processed_passed_reads.map { meta, reads -> tuple(meta, reads, params.downsample_n_reads) }
            ).reads.set { ch_processed_optionally_sampled_reads }
        ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions)
    } else {
        ch_processed_passed_reads.set{ ch_processed_optionally_sampled_reads }
    }
    if ( params.seqtype in ["map-ont", "lr:hq"] && (params.quality_filtering || params.downsample_n_reads)) {
        // MODULE: run NANOPLOT_PROCESSED_READS
        NANOPLOT_PROCESSED_READS(ch_processed_optionally_sampled_reads)
        ch_versions = ch_versions.mix(NANOPLOT_PROCESSED_READS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(
        NANOPLOT_PROCESSED_READS.out.txt.collect{ it[1] }
        )
    }

    // MODULE: run EMU abundance calculation
    EMU_ABUNDANCE(ch_processed_optionally_sampled_reads, ch_emu_db)
    ch_versions = ch_versions.mix(EMU_ABUNDANCE.out.versions)
    if (params.run_krona) {
        // MODULE: Run krona plot
        EMU_REPORT_TO_KRONA(EMU_ABUNDANCE.out.report)
        ch_versions = ch_versions.mix(EMU_REPORT_TO_KRONA.out.versions.first())
        KRONA_KTIMPORTTEXT(EMU_REPORT_TO_KRONA.out.krona_text)
        ch_versions = ch_versions.mix(
            KRONA_KTIMPORTTEXT.out.versions_krona
                .first()
                .map { process, tool, version -> "\"${process}\":\n    ${tool}: ${version}\n" }
                .collectFile(name: 'krona_versions.yml')
        )
    }

    // MODULE: run translate_taxids
    if (params.make_heatmap) {
        TRANSLATE_TAXIDS(EMU_ABUNDANCE.out.assignment_report, PREPARE_EMU_DB.out.taxonomy)
        ch_versions = ch_versions.mix(TRANSLATE_TAXIDS.out.versions)

        // Module: run assignment_heatmap
        ASSIGNMENT_HEATMAP(TRANSLATE_TAXIDS.out.assignment_translated_report)
        ch_versions = ch_versions.mix(ASSIGNMENT_HEATMAP.out.versions)
    }
    // MODULE: run emu combine-outputs
    ch_emu_combine_input_files = channel.empty()
    // Collect all reports into a single list containing the paths
    ch_emu_combine_input_files = EMU_ABUNDANCE.out.report
        .map { it[1] }  // Extract only the file path from the tuple (meta, path)
        .collect()
        .set { collected_files }
    EMU_COMBINE_OUTPUTS(collected_files)
    ch_versions = ch_versions.mix(EMU_COMBINE_OUTPUTS.out.versions)

    // MODULE: run savont on the same reads as emu
    if (params.run_savont && params.seqtype in ["map-ont", "lr:hq"]) {
        SAVONT_ASV(ch_processed_optionally_sampled_reads)
        ch_versions = ch_versions.mix(SAVONT_ASV.out.versions)

        SAVONT_CLASSIFY(SAVONT_ASV.out.asv_dir, ch_emu_db)
        ch_versions = ch_versions.mix(SAVONT_CLASSIFY.out.versions)

        // MODULE: put the emu and savont calls side by side
        COMPARE_CLASSIFIERS(
            EMU_ABUNDANCE.out.report
                .join(SAVONT_CLASSIFY.out.species)
                .join(SAVONT_CLASSIFY.out.mappings)
        )
        ch_versions = ch_versions.mix(COMPARE_CLASSIFIERS.out.versions)
        COMPARE_CLASSIFIERS.out.comparison
            .map { meta, comparison -> comparison }
            .collectFile(
                name: 'classifier_comparison.tsv',
                keepHeader: true,
                skip: 1,
                sort: true,
                storeDir: "${params.outdir}/results"
            )
    }

    // MODULE: run ctrl_comparison
    if (params.ctrl_1) {
        CTRL_COMPARISON(EMU_COMBINE_OUTPUTS.out.combined_report,EMU_COMBINE_OUTPUTS.out.combined_counts_report)
        ch_versions = ch_versions.mix(CTRL_COMPARISON.out.versions)
    }
    // MODULE: Generate PHYLOSEQ object
    if (params.phyloseq) {
        report_ch = EMU_ABUNDANCE.out.report
        all_reports_ch = report_ch
            .map { meta, path -> path }
            .collect()
        COMBINE_REPORTS (
            all_reports_ch
        )

        PHYLOSEQ_OBJECT (
            COMBINE_REPORTS.out.combinedreport,
            PREPARE_EMU_DB.out.taxonomy
        )

        ch_versions = ch_versions.mix(PHYLOSEQ_OBJECT.out.versions)
    }

    // collect tool versions.
    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT REPORTS & MultiQC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // MODULE: MultiQC

    ch_multiqc_config           = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config    = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true) : channel.empty()
    ch_multiqc_logo             = params.multiqc_logo   ? channel.fromPath(params.multiqc_logo, checkIfExists: true) : channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    ch_workflow_summary         = channel.value(paramsSummaryMultiqc(summary_params))
    ch_methods_description      = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )
    multiqc_report = MULTIQC.out.report.toList()

    // MODULE: Generate master.html from reports

    ch_reads
        .map {
            meta, reads -> meta
        }
        .first()
        .set { ch_meta }
    GENERATE_MASTER_HTML(ch_meta, ch_samplesheet)
    ch_versions = ch_versions.mix(GENERATE_MASTER_HTML.out.versions)

    ASSIGNMENT_HEATMAP.out.assignment_heatmap // Use last process so that GENERATE_YAML runs after everything is finished running
        .map {
            meta, assignment_heatmap -> [ meta, file(outdir).toAbsolutePath().toString() ]
        }
        .set { generate_yaml_input }
    GENERATE_YAML(generate_yaml_input)

    emit:
    master_html             = GENERATE_MASTER_HTML.out.html      // channel: [ path(master.html) ]
    versions                = ch_versions                        // channel: [ path(versions.yml) ]
    nanostats_unprocessed   = (params.seqtype in ["map-ont", "lr:hq"]) ? NANOPLOT_UNPROCESSED_READS.out.txt : channel.empty()  // channel: [ path(master.html) ]
    nanostats_processed     = (params.seqtype in ["map-ont", "lr:hq"] && (params.quality_filtering || params.downsample_n_reads)) ? NANOPLOT_PROCESSED_READS.out.txt   : channel.empty()  // channel: [ path(master.html) ]
    multiqc_report          = multiqc_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// An empty gzip file is about 20 bytes
def hasReads(reads) {
    return [reads].flatten().every { it.size() > 20 }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
