process ENTREZDIRECT_EFETCH {
    tag "$sampleid"
    label 'setting_1'
    publishDir { "${params.outdir}/${sampleid}/09_mapping_to_ref" }, mode: 'copy', pattern: '*fasta'

    input:
    tuple val(sampleid), path(ids_to_retrieve)

    output:
    path("*fasta"), optional: true
    tuple val(sampleid), path("${sampleid}_ref_sequences.fasta"), emit: fasta_files, optional: true
    
    script:
    """
    if [ -s "${ids_to_retrieve}" ]; then
      cut -f1 "${ids_to_retrieve}" | while read -r i; do
        efetch -db nucleotide -id "\$i" -format fasta >> "${sampleid}_ref_sequences.fasta"
      done
    fi
    """
}