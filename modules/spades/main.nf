/*
process SPADES {
  tag "${sampleid}"
  label "setting_5"
  containerOptions "${bindOptions}"
  publishDir "${params.outdir}/${sampleid}/06_assembly", mode: 'copy'

  input:
    tuple val(sampleid), path(fastq1), path(fastq2)
  output:
    path("${sampleid}_scaffolds.fasta")
    tuple val(sampleid), path("${sampleid}_scaffolds.fasta"), emit: assembly

  script:
  """
  rnaviralspades.py -1 ${fastq1} \
               -2 ${fastq2} \
               -m 60 -t ${task.cpus} -o ${sampleid}
  cp ${sampleid}/scaffolds.fasta ${sampleid}_scaffolds.fasta
  """
}
*/

process SPADES {
    tag "$meta.id"
    label 'setting_5'

    publishDir "${params.outdir}/${meta.id}/06_assembly", mode: 'copy'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7b/7b7b68c7f8471d9111841dbe594c00a41cdd3b713015c838c4b22705cfbbdfb2/data' :
        'community.wave.seqera.io/library/spades:4.1.0--77799c52e1d1054a' }"

    input:
    tuple val(meta), path(reads)

    output:
    path("*_scaffolds.fasta")
    path("*_contigs.fasta.gz")
    path('*_spades.log')
    path('*_warnings.log'), optional:true
    tuple val(meta.id), path('*_scaffolds.fasta')    , optional:true, emit: assembly
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory.toGiga()
    def readstr = meta.single_end ? "--12 ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    rnaviralspades.py \\
        $args \\
        --threads $task.cpus \\
        --memory $maxmem \\
        ${readstr} \\
        -o spades

    if [ -f spades/scaffolds.fasta ]; then
        mv spades/scaffolds.fasta ${prefix}_scaffolds.fasta
    fi
    if [ -f spades/contigs.fasta ]; then
        mv spades/contigs.fasta ${prefix}_contigs.fasta
        gzip -n ${prefix}_contigs.fasta
    fi
    if [ -f spades/assembly_graph_with_scaffolds.gfa ]; then
        mv spades/assembly_graph_with_scaffolds.gfa ${prefix}_assembly.gfa
        gzip -n ${prefix}_assembly.gfa
    fi
    if [ -f spades/warnings.log ]; then
        mv spades/warnings.log ${prefix}_warnings.log
    fi
    if [ -f spades/spades.log ]; then
        mv spades/spades.log ${prefix}_spades.log
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaviralspades: \$(rnaviralspades.py --version 2>&1 | sed -n 's/^.*SPAdes genome assembler v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_scaffolds.fasta.gz
    echo "" | gzip > ${prefix}_contigs.fasta.gz
    echo "" | gzip > ${prefix}_assembly.gfa.gz
    touch ${prefix}_spades.log
    touch ${prefix}_warnings.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaviralspades: \$(rnaviralspades.py --version 2>&1 | sed -n 's/^.*SPAdes genome assembler v//p')
    END_VERSIONS
    """
}