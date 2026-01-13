/*
//Original process:
//NOTE, the line containerOptions "${bindOptions}" does not work with nfcore modules 
//Returns the following error:
// FATAL:   While checking container encryption: could not open image null: failed to retrieve path for null: lstat null: no such file or directory

process KAIJU {
  tag "${sampleid}"
  label "setting_28"
  containerOptions "${bindOptions}"
  publishDir "${params.outdir}/${sampleid}/05_read_classification", mode: 'copy'

  input:
    tuple val(sampleid), path(fastq1), path(fastq2)
  output:
    file "${sampleid}_kaiju_*name.tsv"
    file "${sampleid}_kaiju_summary*.tsv"
    file "${sampleid}_kaiju.krona"
    tuple val(sampleid), path("${sampleid}_kaiju_summary_viral.tsv"), emit: kaiju_results
    tuple val(sampleid), path ("${sampleid}_kaiju_summary.tsv"), emit: kaiju_results2
    tuple val(sampleid), path("*kaiju.krona"), emit: krona_results


  script:
  """
  c1grep() { grep "\$@" || test \$? = 1; }

  kaiju \
      -z ${task.cpus} \
      -t ${params.kaiju_nodes}  \
      -f ${params.kaiju_dbname} \
      -o ${sampleid}_kaiju.tsv \
      -i ${fastq1} \
      -j ${fastq2} \
      -v

  kaiju-addTaxonNames -t ${params.kaiju_nodes} -n ${params.kaiju_names} -i ${sampleid}_kaiju.tsv -o ${sampleid}_kaiju_name.tsv
  
  kaiju-addTaxonNames \
  -t ${params.kaiju_nodes} \
  -n ${params.kaiju_names} \
  -i ${sampleid}_kaiju.tsv \
  -o ${sampleid}_kaiju_full_lineage_name.tsv \
  -r superkingdom,phylum,class,order,family,genus,species
  
  
  kaiju2table -e -t ${params.kaiju_nodes} -r species -n ${params.kaiju_names} -o ${sampleid}_kaiju_summary.tsv ${sampleid}_kaiju.tsv
  kaiju2krona -t ${params.kaiju_nodes} -n ${params.kaiju_names} -i ${sampleid}_kaiju.tsv -o ${sampleid}_kaiju.krona

  c1grep "taxon_id\\|virus\\|viroid\\|viricota\\|viridae\\|viriform\\|virales\\|virinae\\|viricetes\\|virae\\|viral" ${sampleid}_kaiju_summary.tsv > ${sampleid}_kaiju_summary_viral.tsv
  awk -F'\\t' '\$2>=0.05' ${sampleid}_kaiju_summary_viral.tsv > ${sampleid}_kaiju_summary_viral_filtered.tsv
  """
}

*/
process KAIJU_KAIJU {
    tag "$meta.id"
    label "setting_28"
    

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kaiju:1.10.0--h43eeafb_0':
        'biocontainers/kaiju:1.10.0--h43eeafb_0' }"
    
    input:
    tuple val(meta), path(reads)
    path(db)

    output:
    //tuple val(meta), path('*.tsv'), emit: results
    tuple val(meta.id), path ("${meta.id}_kaiju_summary.tsv"), emit: kaiju_results
    tuple val(meta.id), path("*kaiju.krona"), emit: krona_results
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "-i ${reads}" : "-i ${reads[0]} -j ${reads[1]}"
    """
    dbnodes=`find -L ${db} -name "*nodes.dmp"`
    dbname=`find -L ${db} -name "*kaiju_db_nr_euk.fmi" -not -name "._*"`
    dbnames=`find -L ${db} -name "*names.dmp"`
    

    kaiju \\
        $args \\
        -z ${task.cpus} \\
        -t \$dbnodes\\
        -f \$dbname  \\
        -o ${prefix}_kaiju.tsv \\
        $input

    c1grep() { grep "\$@" || test \$? = 1; }
    kaiju-addTaxonNames -t \$dbnodes -n \$dbnames -i ${prefix}_kaiju.tsv -o ${prefix}_kaiju_name.tsv

    kaiju2table -e -t \$dbnodes -r species -n \$dbnames -o ${prefix}_kaiju_summary.tsv ${prefix}_kaiju.tsv
    kaiju2krona -t \$dbnodes -n \$dbnames -i ${prefix}_kaiju.tsv -o ${prefix}_kaiju.krona

    c1grep "taxon_id\\|virus\\|viroid\\|viricota\\|viridae\\|viriform\\|virales\\|virinae\\|viricetes\\|virae\\|viral" ${prefix}_kaiju_summary.tsv > ${prefix}_kaiju_summary_viral.tsv
    awk -F'\\t' '\$2>=0.05' ${prefix}_kaiju_summary_viral.tsv > ${prefix}_kaiju_summary_viral_filtered.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaiju: \$(echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "-i ${reads}" : "-i ${reads[0]} -j ${reads[1]}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaiju: \$(echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' ))
    END_VERSIONS

    """

}

