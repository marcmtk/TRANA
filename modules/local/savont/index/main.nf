// Build the index SAVONT_CLASSIFY loads, from the database with the exclusions already applied
process SAVONT_INDEX {
    label 'process_low'

    // The fork isn't on bioconda, so there is no conda directive
    // ghcr.io/marcmtk/trana-savont:f828ffd, pinned by digest because the tag could be pushed again
    container 'ghcr.io/marcmtk/trana-savont@sha256:174af1b032ff86c6443c9bfc99d7cae4c9bda9282d633caa547351094e434031'

    input:
    path db

    output:
    path "species_taxid.fasta.mmi" , emit: index
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SAVONT_INDEX needs a container: the fork of savont is not on bioconda. Use a container profile or --run_savont false."
    }
    """
    # savont has no index command, so classify nothing and keep the index it builds. Building
    # it with savont guarantees the settings match, and savont writes it beside the FASTA,
    # hence the private folder: the database itself must not change, or -resume reruns emu.
    mkdir savont_db no_asvs
    ln -s ../$db/.savont_db ../$db/species_taxid.fasta ../$db/taxonomy.tsv savont_db/
    touch no_asvs/final_asvs.fasta

    savont classify \\
        --threads $task.cpus \\
        --input-dir no_asvs \\
        --output-dir no_asvs_classified \\
        --db savont_db

    mv savont_db/species_taxid.fasta.mmi .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        savont: \$(savont --version | sed 's/savont //')-marcmtk-\${SAVONT_COMMIT:0:7}
    END_VERSIONS
    """

    stub:
    """
    touch species_taxid.fasta.mmi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        savont: stub
    END_VERSIONS
    """
}
