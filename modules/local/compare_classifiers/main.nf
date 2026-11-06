// Put the emu and savont calls for one sample side by side
process COMPARE_CLASSIFIERS {
    tag "$meta.id"
    label 'process_single'

    conda 'modules/local/compare_classifiers/env.yaml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/6325e0e1cde97907c30be26d33a170477cca271da36267bf0679965d1b08f0d2/data':
        'community.wave.seqera.io/library/pandas_python:4330fd07d14e9bfb' }"

    input:
    tuple val(meta), path(emu_report), path(savont_species), path(savont_mappings), path(savont_tied), path(emu_check)

    output:
    tuple val(meta), path("${prefix}.classifier_comparison.tsv") , emit: comparison
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def check = emu_check ? "--emu-check ${emu_check}" : ''
    """
    compare_classifiers.py \\
        --sample ${meta.id} \\
        --emu $emu_report \\
        --savont-species $savont_species \\
        --savont-mappings $savont_mappings \\
        --savont-tied $savont_tied \\
        $check \\
        --output ${prefix}.classifier_comparison.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        compare_classifiers.py: \$(compare_classifiers.py --version | sed 's/compare_classifiers.py //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.classifier_comparison.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        compare_classifiers.py: stub
    END_VERSIONS
    """
}
