process RETRIEVE_VIRAL_READS_KRAKEN2 {
    tag "$meta.id"
    label 'setting_11'
    publishDir { "${params.outdir}/${meta.id}/05_read_classification" }, mode: 'copy'

    input:
    tuple val(meta), path(kraken_report), path(kraken_output), path(fastq1), path(fastq2), path(unc_fastq1), path(unc_fastq2)
    output:
    tuple val(meta), path("*_cand_path_R*.fastq.gz"), emit: fastq

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_kraken_reads.py -k ${kraken_output} -r ${kraken_report} \\
                            -t 10239 --include-children \\
                            -s1 ${fastq1} -s2 ${fastq2} \\
                            --fastq-output \\
                            -o ${prefix}_extracted_reads1.fastq -o2 ${prefix}_extracted_reads2.fastq
    gzip ${prefix}_extracted_reads1.fastq
    gzip ${prefix}_extracted_reads2.fastq
    cat ${unc_fastq1} ${prefix}_extracted_reads1.fastq.gz > ${prefix}_cand_path_R1.fastq.gz
    cat ${unc_fastq2} ${prefix}_extracted_reads2.fastq.gz >  ${prefix}_cand_path_R2.fastq.gz
    """
}