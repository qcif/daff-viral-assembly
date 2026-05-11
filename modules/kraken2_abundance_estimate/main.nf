process KRAKEN2_ABUNDANCE_ESTIMATE {
    tag "$meta.id"
    label 'setting_1'

    input:
    tuple val(meta), path(kraken2_report)

    output:
    file("*_kraken_read_abundance.txt")
    tuple val(meta.id), path("*_kraken_read_abundance.txt"), emit: kraken2_results

    when:
    task.ext.when == null || task.ext.when
  
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kraken_lowest_rank.py -i ${kraken2_report} \\
                    -t 3 \\
                    -o ${prefix}_kraken_read_abundance.txt
    """
}