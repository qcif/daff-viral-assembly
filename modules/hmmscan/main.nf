process HMMSCAN {
    tag "${sampleid}"
    label 'setting_20'
    publishDir { "${params.outdir}/${sampleid}/07_annotation" }, mode: 'copy'
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/07/07c4cbd91c4459dc86b13b5cd799cacba96b27d66c276485550d299c7a4c6f8a/data' :
    //   'community.wave.seqera.io/library/hmmer:3.4--cb5d2dd2e85974ca' }"
    containerOptions "--bind ${file(params.hmmer_db).parent}"
    
    input:
    tuple val(sampleid), path(fasta)
    val(hmmer_db)
    
    output:
    file "${sampleid}_orfs.fasta"
    file "${sampleid}_hmmscan*_output.txt"
    tuple val(sampleid), path("${sampleid}_hmmscan_per_target_output.txt"), emit: hmmscan_preds
    tuple val(sampleid), path("${sampleid}_hmmscan_per_domain_output.txt"), emit: hmmscan_domain_preds

    script:
    """
    hmmscan --cpu ${task.cpus} \\
            --domtblout ${sampleid}_hmmscan_per_domain_output.txt \\
            --tblout ${sampleid}_hmmscan_per_target_output.txt \\
            --pfamtblout ${sampleid}_hmmscan_succinct_output.txt \\
            ${hmmer_db} ${fasta} \\
            > ${sampleid}_hmmscan.log 2>&1
    """
}
