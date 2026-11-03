process GENERATE_MASTER_HTML {
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "pyhdfd78af_1") must be EXCLUDED to support installation on different operating systems.
    conda 'modules/local/generate_master_html/env.yaml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nf-core:3.0.2--pyhdfd78af_1':
        'quay.io/biocontainers/nf-core:3.0.2--pyhdfd78af_1' }"

    input:
    val meta
    path csv
    path call_check

    output:
    path output         , emit: html
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    def check = call_check ? "--call-check ${call_check}" : ''
    def prefix = meta.sequencing_run ?: "${params.trace_timestamp}"
    output = "${prefix}_master.html"
    """
    generate_master_html.py --csv ${csv} --html ${params.master_template} --timestamp ${params.trace_timestamp} --output ${output} $args ${check}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_master_html.py: \$(echo \$(generate_master_html.py --version 2>&1) | sed 's/^.*generate_master_html.py //; s/ .*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = meta.sequencing_run ?: "${params.trace_timestamp}"
    output = "${prefix}_master.html"
    """
    touch ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_master_html.py: \$(echo \$(generate_master_html.py --version 2>&1) | sed 's/^.*generate_master_html.py //; s/ .*\$//' )
    END_VERSIONS
    """
}
