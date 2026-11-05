// Put every read of a sample in a bucket: where the pipeline dropped it, what emu made of it and,
// with --structure, how it lies on 16S (full length, truncated, fold-back, repeat, chimera)
process READ_BUCKETS {
    tag "$meta.id"
    label 'process_low'

    // The script needs python, numpy and minimap2, so it runs in emu's own container
    conda 'modules/local/read_buckets/sample/env.yaml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emu:3.6.2--hdfd78af_0':
        'quay.io/biocontainers/emu:3.6.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads, stageAs: 'input/*'), path(dehosted, stageAs: 'dehosted/*'), path(filtered, stageAs: 'filtered/*'), path(sampled, stageAs: 'sampled/*'), path(abundance, stageAs: 'emu/*'), path(assignments, stageAs: 'emu/*'), path(asvs, stageAs: 'savont/*'), path(asv_mappings, stageAs: 'savont/*')
    path db
    path frame

    output:
    tuple val(meta), path("${prefix}.read_buckets.tsv.gz")      , emit: table
    tuple val(meta), path("${prefix}.read_buckets.json")        , emit: summary
    tuple val(meta), path("${prefix}.read_buckets_refs.fasta")  , emit: references
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def optional = [
        dehosted    ? "--dehosted ${dehosted}"       : '',
        filtered    ? "--filtered ${filtered}"       : '',
        sampled     ? "--emu-reads ${sampled}"       : '',
        abundance   ? "--abundance ${abundance}"     : '',
        assignments ? "--assignments ${assignments}" : '',
        asvs && asv_mappings ? "--asvs ${asvs} --asv-mappings ${asv_mappings}" : '',
    ].findAll { it }.join(' ')
    """
    read_buckets.py \\
        --sample ${meta.id} \\
        --prefix ${prefix} \\
        --reads ${reads} \\
        ${optional} \\
        --db ${db} \\
        --frame ${frame} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        read_buckets.py: \$(read_buckets.py --version | sed 's/read_buckets.py //')
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.read_buckets.tsv.gz
    echo '{}' > ${prefix}.read_buckets.json
    touch ${prefix}.read_buckets_refs.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        read_buckets.py: stub
    END_VERSIONS
    """
}
