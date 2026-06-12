process BWA_MEM {
    tag "${sampleid}"
    label 'setting_10'

    input:
    tuple val(sampleid), val(aln_type), path(target), path(fastq1), path(fastq2)

    output:
    tuple val(sampleid), path(target), path("${sampleid}_${aln_type}.sam"), emit: aligned_sam

    script:
    """
    bwa index ${target}
    bwa mem -t ${task.cpus} ${target} $fastq1 $fastq2 > ${sampleid}_${aln_type}.sam 2>> ${sampleid}_${aln_type}.log
    """
}