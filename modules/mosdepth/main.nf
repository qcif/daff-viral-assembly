process MOSDEPTH {
  tag "$sampleid"
  label 'setting_6'

  input:
    tuple val(sampleid), path(consensus), path(bam), path(bai), path(bed)

  output:
    tuple val(sampleid), path("${sampleid}.thresholds.bed"), path("${sampleid}.per-base.bed"), emit: mosdepth_results

  script:
  """
  if [[ ! -s "${consensus}" || ! -s "${bam}" || ! -s "${bai}" || ! -s "${bed}" ]]; then
    touch "${sampleid}.thresholds.bed"
    touch "${sampleid}.per-base.bed"

  else
    mosdepth --by ${bed} --thresholds 1,30  -t ${task.cpus} ${sampleid} ${bam}
    gunzip *.per-base.bed.gz
    gunzip *.thresholds.bed.gz
  fi
  """
}