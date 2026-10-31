// Classify ASVs against the prepared emu database with the marcmtk fork of savont
process SAVONT_CLASSIFY {
    tag "$meta.id"
    label 'process_medium'

    // The fork isn't on bioconda, so there is no conda directive
    container 'ghcr.io/marcmtk/trana-savont:11a4c43'

    input:
    tuple val(meta), path(asv_dir)
    path db

    output:
    tuple val(meta), path("${prefix}.species_abundance.tsv") , emit: species
    tuple val(meta), path("${prefix}.genus_abundance.tsv")   , emit: genus
    tuple val(meta), path("${prefix}.asv_mappings.tsv")      , emit: mappings
    path "${prefix}_savont_classify.log"                     , emit: log
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SAVONT_CLASSIFY needs a container: the fork of savont is not on bioconda. Use a container profile or --run_savont false."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    savont classify \\
        $args \\
        --threads $task.cpus \\
        --input-dir $asv_dir \\
        --output-dir classified \\
        --db $db \\
        > ${prefix}_savont_classify.log 2>&1

    for table in species_abundance genus_abundance asv_mappings; do
        cp classified/\${table}.tsv ${prefix}.\${table}.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        savont: \$(savont --version | sed 's/savont //')-marcmtk-\${SAVONT_COMMIT:0:7}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.species_abundance.tsv ${prefix}.genus_abundance.tsv ${prefix}.asv_mappings.tsv ${prefix}_savont_classify.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        savont: stub
    END_VERSIONS
    """
}
