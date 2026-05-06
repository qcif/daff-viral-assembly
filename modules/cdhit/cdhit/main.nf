process CDHIT_CDHIT {
    tag "${sampleid}"
    label "setting_21"

    input:
    tuple val(sampleid), path(ref)

    output:
    tuple val(sampleid), path("${sampleid}_ref_sequences_clustered.fasta"), emit: clusters

    script:
    """
    cd-hit -i ${ref} -o ${sampleid}_ref_sequences_clustered.fasta -c 0.97 -n 5
    """
}