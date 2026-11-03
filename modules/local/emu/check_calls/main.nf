// Label each emu call confirmed or unconfirmed by the reads that match its references almost exactly
process EMU_CHECK_CALLS {
    tag "$meta.id"
    label 'process_single'

    // The script needs only python and pysam, so it runs in emu's own container
    conda 'modules/local/emu/check_calls/env.yaml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emu:3.6.2--hdfd78af_0':
        'quay.io/biocontainers/emu:3.6.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(report), path(assignments), path(alignments)

    output:
    tuple val(meta), path("${prefix}.emu_call_check.tsv") , emit: check
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    check_emu_calls.py \\
        --sample ${meta.id} \\
        --report $report \\
        --assignments $assignments \\
        --alignments $alignments \\
        --output ${prefix}.emu_call_check.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        check_emu_calls.py: \$(check_emu_calls.py --version | sed 's/check_emu_calls.py //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.emu_call_check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        check_emu_calls.py: stub
    END_VERSIONS
    """
}
