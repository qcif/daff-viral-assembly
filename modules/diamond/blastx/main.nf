process DIAMOND_BLASTX {
    tag "${sampleid}"
    label "setting_30"
    publishDir { "${params.outdir}/${sampleid}/07_annotation" }, mode: 'copy'
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
          ? 'https://depot.galaxyproject.org/singularity/diamond:2.1.24--hf93d47f_0'
          : 'quay.io/biocontainers/diamond:2.1.24--hf93d47f_0'}"

    input:
    tuple val(sampleid), path(viral_fasta), path(other_fasta)
    path(prot_db)
    
    output:
    file "${sampleid}_diamond_matches*.txt"
    tuple val(sampleid), path("${sampleid}_diamond_matches.txt"), emit: diamond_results

    script:
    """
    cat ${viral_fasta} ${other_fasta} > ${sampleid}_combined_contigs.fasta
    diamond blastx --query ${sampleid}_combined_contigs.fasta \\
                  --db ${prot_db} \\
                  --out ${sampleid}_diamond_matches.txt \\
                  --evalue 1e-3 \\
                  --max-target-seqs 1 \\
                  --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle \\
                  --threads ${task.cpus}
    """
}