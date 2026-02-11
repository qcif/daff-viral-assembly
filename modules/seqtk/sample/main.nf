process SEQTK_SAMPLE {
    tag "$meta.id"
    label 'setting_30'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1' :
        'biocontainers/seqtk:1.4--he4a0461_1' }"

    input:
    tuple val(meta), path(reads), path(read_count)
    val(sample_size)
    
    output:
    tuple val(meta), path("*_subsampled.fastq.gz"), emit: reads
    path "versions.yml"                , emit: versions
    when:
    task.ext.when == null || task.ext.when


    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!(args ==~ /.*\ -s\ ?[0-9]+.*/)) {
        args += " -s100"
    }
    if ( !sample_size ) {
        error "SEQTK/SAMPLE must have a sample_size value included"
    }
    """
    # Read count
    READS=\$(cat $read_count)
    THRESHOLD=$sample_size
    if [ "\$READS" -gt "\$THRESHOLD" ]; then
        for f in $reads; do
            gunzip -c "\$f" |
            seqtk \\
                sample \\
                $args \\
                - \\
                $sample_size \\
                | gzip --no-name > "\${f%.fastq.gz}_subsampled.fastq.gz"
        done
    else
        # Just copy or symlink if below threshold
        for f in $reads; do
            cp "\$f" "\${f%.fastq.gz}_subsampled.fastq.gz"
        done
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk 2>&1 | sed -n 's/^Version: //p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.fastq.gz
    """

}