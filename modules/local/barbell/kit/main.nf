// Demultiplex and trim a whole ONT run with the marcmtk fork of barbell
process BARBELL_KIT {
    label 'process_high'

    // The fork isn't on bioconda, so there is no conda directive
    // ghcr.io/marcmtk/trana-barbell:086ab1d, pinned by digest because the tag could be pushed again
    container 'ghcr.io/marcmtk/trana-barbell@sha256:c0ffeb1b368e33dec0323bc8b2f5b142b8b9509fa994f4718316e65ab2e236d0'

    input:
    path fastq_pass
    path barcodes_samplesheet

    output:
    path('fastq_pass_demuxed/*.fastq.gz')     , emit: fastq_files
    path('fastq_pass_demuxed')                , emit: fastq_dir
    path('barbell_read_counts.tsv')           , emit: read_counts
    path('samples_without_reads.txt')         , emit: samples_without_reads
    path('barbell_{annotation,filtered,pattern_per_read}.tsv.gz') , emit: annotation
    path('barbell_log.log')                   , emit: log
    path('versions.yml')                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "BARBELL_KIT needs a container: the fork of barbell is not on bioconda. Use a container profile or --skip_barbell."
    }
    def args = task.ext.args ?: ''
    def sheet = barcodes_samplesheet ? "${barcodes_samplesheet}" : ''
    """
    find -L ${fastq_pass} -name '*.f*q.gz' | sort > input_files.txt
    if [ ! -s input_files.txt ]; then
        echo "No .fastq.gz files found under ${fastq_pass}" >&2
        exit 1
    fi

    barbell kit \\
        $args \\
        --threads $task.cpus \\
        --gzip \\
        --output barbell_out \\
        --input \$(cat input_files.txt) \\
        > barbell_log.log 2>&1

    for table in annotation filtered pattern_per_read; do
        gzip -c barbell_out/\${table}.tsv > barbell_\${table}.tsv.gz
    done

    barbell_bins_to_samples.sh input_files.txt barbell_out fastq_pass_demuxed $sheet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        barbell: \$(barbell --version 2>&1 | grep -o 'barbell [0-9.]*' | sed 's/barbell //')-marcmtk-\${BARBELL_COMMIT:0:7}
    END_VERSIONS
    """

    stub:
    """
    mkdir fastq_pass_demuxed
    echo | gzip > fastq_pass_demuxed/stub.fastq.gz
    printf 'barcode\\tsample_id\\treads\\n' > barbell_read_counts.tsv
    touch samples_without_reads.txt barbell_log.log
    for table in annotation filtered pattern_per_read; do
        echo | gzip > barbell_\${table}.tsv.gz
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        barbell: stub
    END_VERSIONS
    """
}
