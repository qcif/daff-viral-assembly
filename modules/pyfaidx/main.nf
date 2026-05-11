process PYFAIDX {
    tag "$sampleid"
    label 'setting_1'

    input:
    tuple val(sampleid), path(fasta), val(mode)

    output:
    tuple val(sampleid), path("${sampleid}*.bed"), emit: bed

    script:
    def bed_name = mode == 'contigs' ? "${sampleid}_contigs.bed" : "${sampleid}.bed"
    """
    if [[ "${mode}" != "contigs" && "${mode}" != "reference" ]]; then
        echo "Unsupported mode: ${mode}" >&2
        exit 1
    fi

    if [[ ! -s ${fasta} ]]; then
            touch ${bed_name}
    else
            faidx --transform bed ${fasta} > ${bed_name}
    fi
    """
}