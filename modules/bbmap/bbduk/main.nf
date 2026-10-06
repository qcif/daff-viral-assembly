//might want to specify parameter k=31 outside of process in the future
process BBMAP_BBDUK {
    tag "$meta.id"
    label 'setting_11'
    publishDir { "${params.outdir}/${meta.id}/04_cleaned" }, mode: 'copy', pattern: '{*bbduk.log}'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"
    
    input:
    tuple val(meta), path(reads), path(read_count)
    val(mounted_db)
    path(staged_db)
    val(subsample_enabled)
    val(sample_size)

    output:
    //path("${meta.id}_non_rRNA_1.fastq.gz")
    //path("${meta.id}_non_rRNA_2.fastq.gz")
    path("${meta.id}_bbduk.log")
    tuple val(meta), path('*_subsampled*.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    path('*.log')                      , emit: log2
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def trimmed  = meta.single_end ? "out=${prefix}_non_rRNA.fastq.gz" : "out1=${prefix}_non_rRNA_1.fastq.gz out2=${prefix}_non_rRNA_2.fastq.gz"
    def subsampled = meta.single_end ? "out=${prefix}_subsampled.fastq.gz" : "out1=${prefix}_subsampled_1.fastq.gz out2=${prefix}_subsampled_2.fastq.gz"
    def db = mounted_db ?: staged_db
    def contaminants_fa = db ? "ref=${db}" : ''
    if ( subsample_enabled && !sample_size ) {
        error "The bbmap_bbduk process must have a sample_size value included"
    }
    """
    bbduk.sh \\
        -Xmx${task.memory.toGiga()}g \\
        $raw \\
        $trimmed \\
        k=31 \\
        threads=${task.cpus} \\
        $args \\
        $contaminants_fa \\
        &>${prefix}_bbduk.log

    
    if [ "${subsample_enabled}" = "true" ]; then
        # Read count
        READS=\$(cat $read_count)
        THRESHOLD=$sample_size

        #If read counts exceed the threshold, perform subsampling
        if [ "\$READS" -gt "\$THRESHOLD" ]; then
            bbduk.sh \\
            -Xmx${task.memory.toGiga()}g \\
            $trimmed \\
            $subsampled \\
            samplerate=\$(echo "scale=6; \$THRESHOLD / \$READS" | bc) \\
            sampleseed=100
        else
            if [ "$meta.single_end" = true ]; then
                ln -s "${prefix}_non_rRNA.fastq.gz" "${prefix}_subsampled.fastq.gz"
            else
                ln -s "${prefix}_non_rRNA_1.fastq.gz" "${prefix}_subsampled_1.fastq.gz"
                ln -s "${prefix}_non_rRNA_2.fastq.gz" "${prefix}_subsampled_2.fastq.gz"
            fi
        fi
    else
        if [ "$meta.single_end" = true ]; then
            ln -s "${prefix}_non_rRNA.fastq.gz" "${prefix}_subsampled.fastq.gz"
        else
            ln -s "${prefix}_non_rRNA_1.fastq.gz" "${prefix}_subsampled_1.fastq.gz"
            ln -s "${prefix}_non_rRNA_2.fastq.gz" "${prefix}_subsampled_2.fastq.gz"
        fi
    fi
    


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "]")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_command  = meta.single_end ? "echo '' | gzip > ${prefix}.fastq.gz" : "echo '' | gzip > ${prefix}_1.fastq.gz ; echo '' | gzip > ${prefix}_2.fastq.gz"
    """
    touch ${prefix}.bbduk.log
    $output_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "]")
    END_VERSIONS
    """
}