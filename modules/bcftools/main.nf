process BCFTOOLS {
  publishDir { "${params.outdir}/${sampleid}/09_mapping_to_ref" }, mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'

  input:
   tuple val(sampleid), path(ref), path(bam), path(bai)

  output:
    path("${sampleid}_vcf_applied.fasta")
    path("${sampleid}_annotated.vcf.gz")
    tuple val(sampleid), path(ref), path(bam), path(bai), path("${sampleid}_vcf_applied.fasta"), emit: vcf_applied_fasta

  script:
  """
  awk '/^>/ {print; next} {gsub(/[WSMKRYBDHVNwsmskrybdhvn]/, "N"); print}' "${ref}" > "${sampleid}_ref_cleaned.fasta"
  bcftools mpileup -Ou -f ${sampleid}_ref_cleaned.fasta ${bam} | bcftools call -Ou -mv --ploidy=1 | bcftools norm -f ${sampleid}_ref_cleaned.fasta -Oz -o ${sampleid}_raw.vcf.gz
  # -M, --keep-masked-ref           keep sites with masked reference allele (REF=N)
  #-c, --check-ref <e|w|x|s>         check REF alleles and exit (e), warn (w), exclude (x), or set (s) bad sites [e]
  bcftools reheader ${sampleid}_raw.vcf.gz -s <(echo '${sampleid}') \\
  | bcftools filter \\
      -e 'INFO/DP < 20' \\
      -s LOW_DEPTH \\
      --IndelGap 5 \\
      -Oz -o ${sampleid}_annotated.vcf.gz
  
  bcftools index ${sampleid}_annotated.vcf.gz
  # create consensus
  bcftools consensus -f ${sampleid}_ref_cleaned.fasta ${sampleid}_annotated.vcf.gz -o ${sampleid}_vcf_applied.fasta
  """
}