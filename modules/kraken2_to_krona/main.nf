process KRAKEN2_TO_KRONA {
    tag "$meta.id"
    label 'setting_7'
    publishDir { "${params.outdir}/${meta.id}/05_read_classification" }, mode: 'copy'

    input:
    tuple val(meta), path(kraken_report)

    output:
    file("*_kraken_krona.html")

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk '\$1 >= 0.001' ${kraken_report} > filtered_report.txt

    awk -F'\\t' '
    {
        name=\$8
        gsub(/^ +/,"",name)

        indent = match(\$8,/[^ ]/) - 1
        level = int(indent / 2)

        lineage[level] = name

        path = lineage[0]
        for(i=1;i<=level;i++)
            path = path "\\t" lineage[i]

        print \$3 "\\t" path
    }
    ' filtered_report.txt \\
    | ktImportText -o ${prefix}_kraken_krona.html -
    """
}