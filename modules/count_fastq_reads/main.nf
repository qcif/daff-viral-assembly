process COUNT_FASTQ_READS {
  tag "$meta.id"
  label 'setting_1'

  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path(reads), path("*read_count.txt")

  script:
  """
  FWD=\$(ls ${reads} | grep '_1.merged.fastq.gz')
  zgrep -c '^@' "\$FWD" > ${meta.id}_read_count.txt
  """
}