process ORFIPY {
    tag "${sampleid}"
    label "setting_1"
    publishDir { "${params.outdir}/${sampleid}/07_annotation" }, mode: 'copy'

    input:
    tuple val(sampleid), path(viral_fasta), path(other_fasta)
    
    output:
    file "${sampleid}_orfs.fasta"
    tuple val(sampleid), path("${sampleid}_orfs.fasta"), emit: orf_fasta

    script:
    """
    cat ${viral_fasta} ${other_fasta} > ${sampleid}_combined_contigs.fasta
    orfipy ${sampleid}_combined_contigs.fasta \\
      --outdir . \\
      --chunk-size 100000 \\
      --pep ${sampleid}_orfs.fasta \\
      --min 300 \\
      --procs ${task.cpus}
    """
}