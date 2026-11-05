// One report for the run from the per-sample read buckets, barbell's per-read tables and the raw reads
process READ_BUCKETS_REPORT {
    label 'process_single'

    conda 'modules/local/read_buckets/report/env.yaml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emu:3.6.2--hdfd78af_0':
        'quay.io/biocontainers/emu:3.6.2--hdfd78af_0' }"

    input:
    path tables, stageAs: 'samples/*'
    path summaries, stageAs: 'samples/*'
    path barbell_tables, stageAs: 'barbell/*'
    path sheet, stageAs: 'sheet/*'
    path call_check, stageAs: 'check/*'

    output:
    path "read_buckets.html"          , emit: html
    path "read_buckets_summary.tsv"   , emit: summary
    path "read_buckets_calls.tsv"     , emit: calls
    path "*_mqc.tsv"                  , emit: mqc
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def barbell = barbell_tables ? [
        "--barbell-annotation barbell/barbell_annotation.tsv.gz",
        "--barbell-patterns barbell/barbell_pattern_per_read.tsv.gz",
        "--barbell-filtered barbell/barbell_filtered.tsv.gz",
        "--barbell-counts barbell/barbell_read_counts.tsv",
    ].join(' ') : ''
    def optional = [
        barbell,
        sheet      ? "--sheet ${sheet}"             : '',
        call_check ? "--call-check ${call_check}"   : '',
    ].findAll { it }.join(' ')
    """
    read_buckets_report.py \\
        --tables samples/*.read_buckets.tsv.gz \\
        --summaries samples/*.read_buckets.json \\
        ${optional} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        read_buckets_report.py: \$(read_buckets_report.py --version | sed 's/read_buckets_report.py //')
    END_VERSIONS
    """

    stub:
    """
    touch read_buckets.html read_buckets_summary.tsv read_buckets_calls.tsv read_buckets_fate_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        read_buckets_report.py: stub
    END_VERSIONS
    """
}
