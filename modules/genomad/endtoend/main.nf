process GENOMAD_ENDTOEND {
    tag "${sampleid}"
    label "setting_27"
    publishDir { "${params.outdir}/${sampleid}/07_annotation/genomad" }, mode: 'copy'

    input:
    tuple val(sampleid), path(viral_fasta), path(other_fasta)
    path(genomad_db)

    output:
    file "*_summary/*_virus.fna"
    file "*_summary/*_virus_summary.tsv"
    file "*_summary/*_virus_genes.tsv"
    file "*_summary/*_virus_proteins.faa"
    file "*_find_proviruses/*_provirus.tsv"
    file "*_find_proviruses/*_provirus_taxonomy.tsv"
    file "*_find_proviruses/*_provirus.fna"
    file "*_find_proviruses/*_provirus_genes.tsv"
    file "*_find_proviruses/*_provirus_proteins.faa"
    file "*_aggregated_classification/*aggregated_classification.tsv"
    file "*_annotate/*_taxonomy.tsv"
    file "${sampleid}_genomad.log"
    tuple val(sampleid), path("${sampleid}_combined_contigs_virus_summary.tsv"), emit: virus_preds
      
    script:
    """
    cat ${viral_fasta} ${other_fasta} > ${sampleid}_combined_contigs.fasta
    genomad \\
      end-to-end \\
      ${sampleid}_combined_contigs.fasta \\
      ./ \\
      ${genomad_db} \\
      --threads ${task.cpus} \\
      --min-score 0.7 \\
      --splits 1 \\
      > ${sampleid}_genomad.log 2>&1
      cp ${sampleid}_combined_contigs_summary/${sampleid}_combined_contigs_virus_summary.tsv .
    """
}