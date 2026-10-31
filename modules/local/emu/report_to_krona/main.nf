// Turn an emu abundance report into the text format of ktImportText
process EMU_REPORT_TO_KRONA {
    tag "$meta.id"
    label 'process_single'

    conda 'modules/local/emu/report_to_krona/env.yaml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/6325e0e1cde97907c30be26d33a170477cca271da36267bf0679965d1b08f0d2/data':
        'community.wave.seqera.io/library/pandas_python:4330fd07d14e9bfb' }"

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path("${prefix}.krona.txt") , emit: krona_text
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    emu_report_to_krona.py $report ${prefix}.krona.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emu_report_to_krona.py: \$(emu_report_to_krona.py --version | sed 's/emu_report_to_krona.py //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.krona.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emu_report_to_krona.py: stub
    END_VERSIONS
    """
}
