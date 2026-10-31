// Unpack the emu database into the one form that both emu and savont can read
process PREPARE_EMU_DB {
    label 'process_single'

    conda 'modules/local/prepare_emu_db/env.yaml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/6325e0e1cde97907c30be26d33a170477cca271da36267bf0679965d1b08f0d2/data':
        'community.wave.seqera.io/library/pandas_python:4330fd07d14e9bfb' }"

    input:
    path db
    path exclude_ids

    output:
    path "emu_db"                , emit: db
    path "emu_db/taxonomy.tsv"   , emit: taxonomy
    path "emu_db_summary.tsv"    , emit: summary
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def exclude = exclude_ids ? "--exclude-ids ${exclude_ids}" : ''
    """
    prepare_emu_db.py \\
        --db $db \\
        --outdir emu_db \\
        --summary emu_db_summary.tsv \\
        $exclude

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prepare_emu_db.py: \$(prepare_emu_db.py --version | sed 's/prepare_emu_db.py //')
    END_VERSIONS
    """

    stub:
    """
    mkdir emu_db
    touch emu_db/species_taxid.fasta emu_db/taxonomy.tsv emu_db_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prepare_emu_db.py: stub
    END_VERSIONS
    """
}
