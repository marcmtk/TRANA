// Call ASVs with the marcmtk fork of savont
process SAVONT_ASV {
    tag "$meta.id"
    label 'process_medium'

    // The fork isn't on bioconda, so there is no conda directive
    container 'ghcr.io/marcmtk/trana-savont:11a4c43'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}_asv")                 , emit: asv_dir
    tuple val(meta), path("${prefix}.final_asvs.fasta")    , emit: fasta
    tuple val(meta), path("${prefix}.feature-table.tsv")   , emit: feature_table
    path "${prefix}_savont_asv.log"                        , emit: log
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SAVONT_ASV needs a container: the fork of savont is not on bioconda. Use a container profile or --run_savont false."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    savont asv \\
        $args \\
        --threads $task.cpus \\
        --output-dir ${prefix}_asv \\
        $reads \\
        > ${prefix}_savont_asv.log 2>&1

    # Only clustering intermediates, and they can get large
    rm -rf ${prefix}_asv/temp

    cp ${prefix}_asv/final_asvs.fasta ${prefix}.final_asvs.fasta
    cp ${prefix}_asv/feature-table.tsv ${prefix}.feature-table.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        savont: \$(savont --version | sed 's/savont //')-marcmtk-\${SAVONT_COMMIT:0:7}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_asv
    touch ${prefix}.final_asvs.fasta ${prefix}.feature-table.tsv ${prefix}_savont_asv.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        savont: stub
    END_VERSIONS
    """
}
