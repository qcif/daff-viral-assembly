process MOSDEPTH {
  tag "$sampleid"
  label "setting_2"

  input:
    tuple val(sampleid), path(consensus), path(bam), path(bai), path(bed)

  output:
    tuple val(sampleid), path("${sampleid}.thresholds.bed"), emit: mosdepth_results

  script:
  """
  if [[ ! -s ${consensus} ]]; then
    touch ${sampleid}.thresholds.bed

  else
    mosdepth --by ${bed} --thresholds 30 -t ${task.cpus} ${sampleid} ${bam}
    gunzip *.per-base.bed.gz
    gunzip *.thresholds.bed.gz
  fi
  """
}