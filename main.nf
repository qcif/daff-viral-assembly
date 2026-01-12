#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { fromSamplesheet } from 'plugin/nf-validation'


def helpMessage () {
    log.info """
    daff-viral-assembly
    Marie-Emilie Gauthier
    Cameron Hyde

    Usage:
    Run the command
    nextflow run main.nf -profile singularity -params-file {params.yml}

    Required arguments:
      --analysis_mode                 clustering, map2ref
                                      Default: 'clustering' [required]
      --analyst_name                  Name of the analyst
                                      Default: null [required]
      --facility                      Name of the facility where the analyst is performing the analysis
                                      Default: null [required]

    Optional arguments:
      --help                          Will print this usage document
      -resume                         Resume a failed run
      --outdir                        Path to save the output file
                                      'results'
      --samplesheet '[path/to/file]'  Path to the csv file that contains the list of
                                      samples to be analysed by this pipeline.
                                      Default:  'index.csv'. Needs to have .csv suffix


    Contents of index.csv:
      sampleid,fastq_path,target_organisms,target_gene,target_size,fwd_primer,rev_primer
      VE24-1279_COI,/work/tests/mtdt_data/barcode01_VE24-1279_COI/*fastq.gz,drosophilidae,COI,711,GGTCAACAAATCATAAAGATATTGG,ATTTTTTGGTCACCCTGAAGTTTA

      #### Blast options ####
      --blast_threads                 Number of threads for megablast
                                      Default: '2'
      --blastn_db                     Path to blast database [required if not performing qc_only or preprocessing_only]
                                      Default: ''
      --taxdump                       Path to taxonomykit database directory [required if not performing qc_only or preprocessing_only]
                                      Default: ''

      #### Mapping back to ref options ####
      --mapping_back_to_ref           Mapped back to reference blast match
                                      Default: 'true'

    """.stripIndent()
}
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}
if (params.blastn_db != null) {
    blastn_db_name = file(params.blastn_db).name
    blastn_db_dir = file(params.blastn_db).parent
}

//if (params.sortmerna_ref != null) {
//    sortmerna_ref_name = file(params.sortmerna_ref).name
//    sortmerna_ref_dir = file(params.sortmerna_ref).parent
//}

if (params.kraken2_db != null) {
    krkdb_dir = file(params.kraken2_db).parent
}
//if (params.kaiju_nodes != null & params.kaiju_dbname != null & params.kaiju_names != null) {
//    kaiju_dbs_dir = file(params.kaiju_nodes).parent
//}
/*
if (params.kaiju_db_path != null ) {
    kaiju_dbs_dir = file(params.kaiju_nodes).parent
}
*/
//if (params.diamond_db != null) {
//    diamond_db_dir = file(params.diamond_db).parent
//}

if (params.genomad_db != null) {
    genomad_db_dir = file(params.genomad_db).parent
}
if (params.hmmer_db != null) {
    hmmer_db_dir = file(params.hmmer_db).parent
}
//if (params.taxdump != null) {
//    taxdump_dir = file(params.taxdump).parent
//}

//if (params.reference != null) {
//    reference_name = file(params.reference).name
 //   reference_dir = file(params.reference).parent
//}
//if (params.host_fasta != null) {
//   host_fasta_dir = file(params.host_fasta).parent
//}

def isNonEmptyFile(file) {
    return file.exists() && file.size() > 0
}

switch (workflow.containerEngine) {
  case "singularity":
    bindbuild = "";
    if (params.blastn_db != null) {
      bindbuild = (bindbuild + "-B ${blastn_db_dir} ")
    }
    if (params.taxdump != null) {
      bindbuild = (bindbuild + "-B ${params.taxdump} ")
    }

    if (params.kaiju_db_path != null) {
      bindbuild = (bindbuild + "-B ${params.kaiju_db_path} ")
    }
//    if (params.diamond_db != null) {
//      bindbuild = (bindbuild + "-B ${diamond_db_dir} ")
//    }

//    if (params.sortmerna_ref != null) {
//      bindbuild = (bindbuild + "-B ${sortmerna_ref_dir} ")
//    }

    if (params.kraken2_db != null) {
      bindbuild = (bindbuild + "-B ${krkdb_dir} ")
    }
    
    //if (params.kaiju_nodes != null & params.kaiju_dbname != null & params.kaiju_names != null) {
    //  bindbuild = (bindbuild + "-B ${kaiju_dbs_dir} ")
    //}
    
    if (params.genomad_db != null ) {
      bindbuild = (bindbuild + "-B ${genomad_db_dir} ")
    }
    if (params.hmmer_db != null ) {
      bindbuild = (bindbuild + "-B ${hmmer_db_dir} ")
    }

//    if (params.reference != null) {
//      bindbuild = (bindbuild + "-B ${reference_dir} ")
//    }
//   if (params.host_fasta != null) {
//      bindbuild = (bindbuild + "-B ${host_fasta_dir} ")
//    }
    bindOptions = bindbuild;
    break;
  default:
    bindOptions = "";
}

process BLASTN {
  tag "${sampleid}"
  containerOptions "${bindOptions}"
  label "setting_20"

  input:
    tuple val(sampleid), path(assembly)
  output:
    tuple val(sampleid), path("${sampleid}*_blastn.bls"), emit: blast_results

  script:
  def blastoutput = assembly.getBaseName() + "_blastn.bls"
    """
    blastn -query ${assembly} \
      -db ${params.blastn_db} \
      -out ${blastoutput} \
      -evalue 1e-3 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
      -max_target_seqs 5
    """
}

process COVSTATS {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/09_mapping_to_ref", mode: 'copy'

  input:
    tuple val(sampleid), path(bed), path(blast_results), path(bbsplit_stats), path(consensus), path(coverage), path(mapping_q), path(ref)
  output:
    path("*reference_with_cov_stats.txt")
    tuple val(sampleid), path("*reference_with_cov_stats.txt"), emit: detections_summary
    tuple val(sampleid), path(consensus), path("*reference_with_cov_stats.txt"), emit: detections_summary3
    path("*reference_with_cov_stats.txt"), emit: detections_summary2


  script:
    """
    derive_coverage_stats.py --sample ${sampleid} --blastn_results ${blast_results} --bbsplit_stats ${bbsplit_stats} --coverage ${coverage} --bed ${bed} --reference ${ref} --consensus ${consensus} --mapping_quality ${mapping_q}
    """
}

process CONTIG_COVSTATS {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/08_mapping_to_contigs", mode: 'copy'

  input:
    tuple val(sampleid), path(bed), path(blast_results), path(bbsplit_stats), path(coverage), path(mapping_q)
  output:
    path("*_with_cov_stats.txt")
    tuple val(sampleid), path("*_with_cov_stats.txt"), emit: detections_summary
//    tuple val(sampleid), path(consensus), path("*reference_with_cov_stats.txt"), emit: detections_summary3
//    path("*reference_with_cov_stats.txt"), emit: detections_summary2


  script:
    """
    derive_contig_coverage_stats.py --sample ${sampleid} --blastn_results ${blast_results} --bbsplit_stats ${bbsplit_stats} --coverage ${coverage} --bed ${bed} --mapping_quality ${mapping_q}
    """
}


process EXTRACT_BLAST_HITS {
  tag "${sampleid}"
  label "setting_1"
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(blast_results), val(target_organism), val(target_gene), val(target_size)

  output:
    tuple val(sampleid), path("${sampleid}*_megablast_top_hits_tmp.txt"), emit: topblast
    tuple val(sampleid), path("${sampleid}_reference_match.fasta"), emit: reference_fasta_files
    tuple val(sampleid), path("${sampleid}_final_polished_consensus_match.fasta"), emit: consensus_fasta_files

  script:
    target_organism_str = (target_organism instanceof List)
    ? "\"${target_organism.join('|')}\""
    : "\"${target_organism}\""
    """
    if [[ \$(wc -l < *_megablast_top_10_hits.txt) -ge 2 ]]
      then
        select_top_blast_hit.py --sample_name ${sampleid} --blastn_results ${sampleid}*_top_10_hits.txt --target_organism ${target_organism_str} --taxonkit_database_dir ${params.taxdump}

        # extract segment of consensus sequence that align to reference
        awk  -F  '\\t' 'NR>1 { printf ">%s\\n%s\\n",\$2,\$23 }' ${sampleid}*_top_hits_tmp.txt | sed 's/-//g' > ${sampleid}_final_polished_consensus_match.fasta

        # extract segment of reference that align to consensus sequence
        awk  -F  '\\t' 'NR>1 { printf ">%s_%s\\n%s\\n",\$2,\$4,\$24 }' ${sampleid}*_top_hits_tmp.txt | sed 's/-//g' > ${sampleid}_reference_match.fasta
    else
        echo "No hits found for ${sampleid} in the blast results. Skipping the extraction of consensus and reference fasta files." >&2
        touch ${sampleid}_final_polished_consensus_match.fasta
        touch ${sampleid}_reference_match.fasta
        touch ${sampleid}_megablast_top_hits_tmp.txt
    fi
    """
}

process FASTA2TABLE {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/07_annotation", mode: 'copy'

  input:
    tuple val(sampleid), path(tophits), path(fasta)
  output:
    file("${sampleid}_megablast_top_viral_hits_with_contigs.txt")
    file("${sampleid}_megablast_top_viral_hits_filtered_with_contigs.txt")
    tuple val(sampleid), file("${sampleid}_ref_ids_to_retrieve.txt"), emit: ref_ids
    tuple val(sampleid), file("${sampleid}_contig_ids_to_retrieve.txt"), emit: contig_ids
    tuple val(sampleid), file("${sampleid}_megablast_top_viral_hits_filtered_with_contigs.txt"), emit: blast_results
    tuple val(sampleid), file("${sampleid}_megablast_top_viral_hits_with_contigs.txt"), emit: blast_results2


  script:
    """
    fasta2table.py --fasta ${fasta} --sample ${sampleid} --tophits ${tophits}
    cut -f2 ${sampleid}_megablast_top_viral_hits_filtered_with_contigs.txt | sed '1d' | sed 's/ //g' | sort | uniq > ${sampleid}_contig_ids_to_retrieve.txt
    """
}

process FASTA2TABLE2 {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/09_mapping_to_ref", mode: 'copy'

  input:
    tuple val(sampleid), path(fasta), path(stats)
  output:
    file("${sampleid}_reference_with_cov_stats_final.txt")
    tuple val(sampleid), file("${sampleid}_reference_with_cov_stats_final.txt"), emit: detections_summary_final

  script:
    """
    fasta2table.py --fasta ${fasta} --sample ${sampleid} --tophits ${stats}
    """
}

process MOSDEPTH {
  tag "$sampleid"
  label "setting_3"

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


process MOSDEPTH_CONTIGS {
  tag "$sampleid"
  label "setting_3"

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



process PYFAIDX {
  tag "$sampleid"
  label "setting_3"

  input:
    tuple val(sampleid), path(fasta)

  output:
    tuple val(sampleid), path("${sampleid}.bed"), emit: bed

  script:
    """
    if [[ ! -s ${fasta} ]]; then
      touch ${sampleid}.bed
    else
      faidx --transform bed ${fasta} > ${sampleid}.bed
    fi
    """
}

process PYFAIDX_CONTIGS {
  tag "$sampleid"
  label "setting_3"

  input:
    tuple val(sampleid), path(fasta)

  output:
    tuple val(sampleid), path("${sampleid}_contigs.bed"), emit: bed

  script:
    """
    if [[ ! -s ${fasta} ]]; then
      touch ${sampleid}.bed
    else
      faidx --transform bed ${fasta} > ${sampleid}_contigs.bed
    fi
    """
}

process QCREPORT {
  publishDir "${params.outdir}/02_qc_report", mode: 'copy', overwrite: true
  containerOptions "${bindOptions}"

  input:
    path multiqc_files

  output:
    path("run_qc_report_*txt")
    path("run_qc_report_*html")
    path("run_qc_report_*html"), emit: qc_report_html
    path("run_qc_report_*txt"), emit: qc_report_txt

  script:
    """
    seq_run_qc_report.py
    """
}

process TIMESTAMP_START {
  publishDir "${params.outdir}/01_pipeline_logs", mode: 'copy', overwrite: true
  cache false
  output:
  path "*nextflow_start_timestamp.txt"
  path("*nextflow_start_timestamp.txt"), emit: timestamp

  script:
    """
    START_TIMESTAMP=\$(date "+%Y%m%d%H%M%S")
    echo "\$START_TIMESTAMP" > "\${START_TIMESTAMP}_nextflow_start_timestamp.txt"
    """
}

process HTML_REPORT {
  publishDir "${params.outdir}/${sampleid}/11_report", mode: 'copy', overwrite: true
  containerOptions "${bindOptions}"
  label 'setting_3'

  input:
    tuple val(sampleid), path(raw_fastqc), path(filtered_fastqc), path (fastp), path(fasta), path(summary_known_viruses), path(kaiju_summary), path(kraken_summary), path(detections_summary), path(mapping_summary), path(consensus), path(bam), path(bai), path(novel_virus_summary),
    path(timestamp),
    path(qcreport_html),
    path(qcreport_txt),
    path(configyaml),
    path(samplesheet)

  output:
    path("*"), optional: true
    path("run_qc_report.html"), optional: true

  script:
  //analyst_name = params.analyst_name.replaceAll(/ /, '_')
  //facility = params.facility.replaceAll(/ /, '_')
  analyst_name = "Maely"
  facility = "QUT"
    """
    cp ${qcreport_html} run_qc_report.html
    cp ${params.tool_versions} versions.yml
    cp ${params.default_params} default_params.yml

    build_report.py --samplesheet ${samplesheet} --result_dir . --params_file ${configyaml} --analyst ${analyst_name} --facility ${facility} --versions versions.yml --default_params_file default_params.yml
    #build_report.py --samplesheet ${samplesheet} --result_dir . --params_file ${configyaml} --versions versions.yml --default_params_file default_params.yml
    
    """
}

process SEQTK {
  tag "${sampleid}"
  label "setting_2"

  input:
    tuple val(sampleid), path(assembly)
  output:
    path("${sampleid}_scaffolds_filt.fasta")
    tuple val(sampleid), path("${sampleid}_scaffolds_filt.fasta"), emit: filt_fasta

  script:
  """
  seqtk seq -L 150 ${assembly} > ${sampleid}_scaffolds_filt.fasta
  
  """
}
/*
process SORTMERNA {
  tag "${sampleid}"
  label "setting_21"
  containerOptions "${bindOptions}"
  publishDir "${params.outdir}/${sampleid}/04_cleaned", mode: 'copy'

  input:
    tuple val(sampleid), path(fastq1), path(fastq2), path(ref), path(index)
  output:
    path("${sampleid}_rRNA_reads.log")
    path("${sampleid}_non_rRNA_reads_fwd.fq.gz")
    path("${sampleid}_non_rRNA_reads_rev.fq.gz")
    tuple val(sampleid), path("${sampleid}_non_rRNA_reads_fwd.fq.gz"), path("${sampleid}_non_rRNA_reads_rev.fq.gz"), emit: fastp_filtered_fq

  script:
  """
  sortmerna \
    --ref ${ref} \
    --reads ${fastq1} \
    --reads ${fastq2} \
    --threads ${task.cpus} --workdir . --aligned ${sampleid}_rRNA_reads \
    --fastx --other ${sampleid}_non_rRNA_reads --paired_in --out2 --num_alignments 1 -v --index 0
  """
}
*/
process SPADES {
  tag "${sampleid}"
  label "setting_5"
  containerOptions "${bindOptions}"
  publishDir "${params.outdir}/${sampleid}/06_assembly", mode: 'copy'

  input:
    tuple val(sampleid), path(fastq1), path(fastq2)
  output:
    path("${sampleid}_scaffolds.fasta")
    tuple val(sampleid), path("${sampleid}_scaffolds.fasta"), emit: assembly

  script:
  """
  rnaviralspades.py -1 ${fastq1} \
               -2 ${fastq2} \
               -m 60 -t ${task.cpus} -o ${sampleid}
  cp ${sampleid}/scaffolds.fasta ${sampleid}_scaffolds.fasta
  """
}

process EXTRACT_VIRAL_CLASSIFICATION_HITS {
  tag "${sampleid}"
  label "setting_2"
  publishDir "${params.outdir}/${sampleid}/07_annotation", mode: 'copy'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(blast_results) 

  output:
    //tuple val(sampleid), path("${sampleid}_megablast_top_viral_hits_filtered.txt"), emit: viral_blast_results
    tuple val(sampleid), path("${sampleid}_megablast_top_viral_hits.txt"), emit: viral_blast_results
    //path("${sampleid}_megablast_top_viral_hits.txt")
    //path("${sampleid}_megablast_top_viral_hits_filtered.txt")
    path("${sampleid}_blastn.txt")

  script:
  """
  cat ${blast_results} > ${sampleid}_blastn.txt
  filter_blast.py --blastn_results ${sampleid}_blastn.txt --sample_name ${sampleid} --taxonkit_database_dir ${params.taxdump} --filter ${params.filter_terms}
  """
}

process EXTRACT_VIRAL_BLAST_HITS {
  tag "${sampleid}"
  label "setting_2"
  publishDir "${params.outdir}/${sampleid}/07_annotation", mode: 'copy'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(blast_results) 

  output:
    //tuple val(sampleid), path("${sampleid}_megablast_top_viral_hits_filtered.txt"), emit: viral_blast_results
    tuple val(sampleid), path("${sampleid}_megablast_top_viral_hits.txt"), emit: viral_blast_results
    //path("${sampleid}_megablast_top_viral_hits.txt")
    //path("${sampleid}_megablast_top_viral_hits_filtered.txt")
    path("${sampleid}_blastn.txt")

  script:
  """
  cat ${blast_results} > ${sampleid}_blastn.txt
  filter_blast.py --blastn_results ${sampleid}_blastn.txt --sample_name ${sampleid} --taxonkit_database_dir ${params.taxdump} --filter ${params.filter_terms}
  """
}

process SUMMARISE_READ_CLASSIFICATION {
  tag "${sampleid}"
  label "setting_2"
  publishDir "${params.outdir}/${sampleid}/05_read_classification", mode: 'copy'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(kaiju_results), path(bracken_results), path(stats)

  output:
    //tuple val(sampleid), path("${sampleid}_megablast_top_viral_hits_filtered.txt"), emit: viral_blast_results
    path("${sampleid}_kaiju_summary.txt")
    path("${sampleid}_kraken_summary.txt")
    tuple val(sampleid), path("${sampleid}_kaiju_summary.txt"), emit: kaiju_summary
    tuple val(sampleid), path("${sampleid}_kraken_summary.txt"), emit: kraken_summary
    //path("${sampleid}_megablast_top_viral_hits.txt")
    //path("${sampleid}_megablast_top_viral_hits_filtered.txt")

  script:
  """
  filter_classification_results.py --kaiju ${kaiju_results} --sample_name ${sampleid} --bracken ${bracken_results} --taxonkit_database_dir ${params.taxdump} --stats ${stats} --filter ${params.filter_terms}
  """
}

process EXTRACT_REF_FASTA {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/09_mapping_to_ref", mode: 'copy', pattern: '*fasta'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(ids_to_retrieve)

  output:
    path("*fasta"), optional: true
    tuple val(sampleid), path("${sampleid}_ref_sequences.fasta"), emit: fasta_files, optional: true
  
  script:
    """
    #cut -f1,4 ${ids_to_retrieve} | sed '1d' | sed 's/ /_/g' | sort | uniq > ids_to_retrieve.txt
    #if [ -s ${ids_to_retrieve} ]
    #  then
    #    for i in `cut -f2  ids_to_retrieve.txt`; do j=`grep \${i} ids_to_retrieve.txt | cut -f1`; efetch -db nucleotide  -id \${i} -format fasta > ${sampleid}_\${i}.fasta ; done
    #    for i in `cut -f1  ${ids_to_retrieve}`; do efetch -db nucleotide  -id \${i} -format fasta > ${sampleid}_\${i}.fasta ; done
    #    cat ${sampleid}_*.fasta > ${sampleid}_ref_sequences.fasta
    #fi

    if [ -s "${ids_to_retrieve}" ]; then
    
       cut -f1 "${ids_to_retrieve}" | while read -r i; do
          efetch -db nucleotide -id "\$i" -format fasta >> "${sampleid}_ref_sequences.fasta"
      done
    fi
    """
}

process CLUSTER {
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

process MAPPING_BACK_TO_REF {
  tag "${sampleid}"
  label "setting_21"

  input:
    tuple val(sampleid), path(ref), path(fastq1), path(fastq2)

  output:
    tuple val(sampleid), path(ref), file("${sampleid}_ref_aln.sam"), emit: aligned_sam

  script:
  """
  bowtie2-build $ref $ref
  bowtie2 --threads ${task.cpus} --very-sensitive-local -k 100 --score-min L,20,1.0  -x $ref \
          -1 $fastq1 -2 $fastq2 -S ${sampleid}_ref_aln.sam 2>> ${sampleid}_mapping.log
  """
}

process MAPPING_BACK_TO_CONTIGS {
  tag "${sampleid}"
  label "setting_21"

  input:
    tuple val(sampleid), path(contigs), path(fastq1), path(fastq2)

  output:
    tuple val(sampleid), path(contigs), file("${sampleid}_contig_aln.sam"), emit: contig_aligned_sam

  script:
  """
  bowtie2-build ${contigs} ${contigs}
  bowtie2 --threads ${task.cpus} --very-sensitive-local -k 100 --score-min L,20,1.0  -x ${contigs} \
          -1 $fastq1 -2 $fastq2 -S ${sampleid}_contig_aln.sam 2>> ${sampleid}_mapping.log
  """
}

process SAMTOOLS2 {
  publishDir "${params.outdir}/${sampleid}/09_mapping_to_ref", mode: 'copy'
  tag "${sampleid}"
  label 'setting_21'

  input:
    tuple val(sampleid), path(ref), path(sam)

  output:
    path "${sampleid}_ref_aln.sorted.bam"
    path "${sampleid}_ref_aln.sorted.bam.bai"
    path "${sampleid}_samtools_consensus_from_ref.fasta"
    path "${sampleid}_ref_coverage.txt"
    tuple val(sampleid), path(ref), path("${sampleid}_ref_aln.sorted.bam"), path("${sampleid}_ref_aln.sorted.bam.bai"), emit: sorted_bam
    tuple val(sampleid), path("${sampleid}_ref_coverage.txt"), emit: coverage
    tuple val(sampleid), path("${sampleid}_ref_mapq.txt"), emit: mapping_quality

  script:
    """
    #filter out unmapped reads -F 4
    samtools view -@ ${task.cpus} -Sb -F 4 ${sam} | samtools sort -@ ${task.cpus} -o ${sampleid}_ref_aln.sorted.bam
    samtools index ${sampleid}_ref_aln.sorted.bam
    samtools coverage ${sampleid}_ref_aln.sorted.bam  > ${sampleid}_ref_coverage.txt
    samtools view -@ ${task.cpus} ${sampleid}_ref_aln.sorted.bam | awk '{mapq[\$3]+=\$5; count[\$3]++} END {for (chr in mapq) printf "%s\\t%.2f\\n", chr, mapq[chr]/count[chr]}' > ${sampleid}_ref_mapq.txt
    samtools consensus -@ ${task.cpus} -f fastq -a -A ${sampleid}_ref_aln.sorted.bam -o ${sampleid}_consensus.fastq 
    samtools consensus -@ ${task.cpus} -f fasta -a -A ${sampleid}_ref_aln.sorted.bam -o ${sampleid}_samtools_consensus_from_ref.fasta
    """
}

process SAMTOOLS_CONTIGS {
  publishDir "${params.outdir}/${sampleid}/08_mapping_to_contigs", mode: 'copy'
  tag "${sampleid}"
  label 'setting_21'

  input:
    tuple val(sampleid), path(ref), path(sam)

  output:
    path "${sampleid}_contig_aln.sorted.bam"
    path "${sampleid}_contig_aln.sorted.bam.bai"
    path "${sampleid}_contig_coverage.txt"
    tuple val(sampleid), path(ref), path("${sampleid}_contig_aln.sorted.bam"), path("${sampleid}_contig_aln.sorted.bam.bai"), emit: sorted_bam
    tuple val(sampleid), path("${sampleid}_contig_coverage.txt"), emit: coverage
    tuple val(sampleid), path("${sampleid}_contig_mapq.txt"), emit: mapping_quality

  script:
    """
    #filter out unmapped reads -F 4
    samtools view -@ ${task.cpus} -Sb -F 4 ${sam} | samtools sort -@ ${task.cpus} -o ${sampleid}_contig_aln.sorted.bam
    samtools index ${sampleid}_contig_aln.sorted.bam
    samtools coverage ${sampleid}_contig_aln.sorted.bam  > ${sampleid}_contig_coverage.txt
    samtools view -@ ${task.cpus} ${sampleid}_contig_aln.sorted.bam | awk '{mapq[\$3]+=\$5; count[\$3]++} END {for (chr in mapq) printf "%s\\t%.2f\\n", chr, mapq[chr]/count[chr]}' > ${sampleid}_contig_mapq.txt
    """
}

process BCFTOOLS {
  publishDir "${params.outdir}/${sampleid}/09_mapping_to_ref", mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'
  containerOptions "${bindOptions}"

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
		bcftools reheader ${sampleid}_raw.vcf.gz -s <(echo '${sampleid}') \
    | bcftools filter \
        -e 'INFO/DP < 20' \
        -s LOW_DEPTH \
        --IndelGap 5 \
        -Oz -o ${sampleid}_annotated.vcf.gz
    
    bcftools index ${sampleid}_annotated.vcf.gz
		# create consensus
    bcftools consensus -f ${sampleid}_ref_cleaned.fasta ${sampleid}_annotated.vcf.gz -o ${sampleid}_vcf_applied.fasta \
    """
}

process BEDTOOLS {
  tag "${sampleid}"
  label 'setting_3'
  containerOptions "${bindOptions}"
  publishDir "${params.outdir}/${sampleid}/09_mapping_to_ref", mode: 'copy'

  input:
   tuple val(sampleid), path(ref), path(bam), path(bai), path(vcf_applied_fasta)

  output:
    path("${sampleid}_bcftools_masked_consensus.fasta")
     tuple val(sampleid), path("${sampleid}_bcftools_masked_consensus.fasta"), emit: bcftools_masked_consensus_fasta

  script:
    """
		bedtools genomecov -ibam ${bam} -bga > ${sampleid}_genomecov.bed
    awk '\$4==0 {print}' ${sampleid}_genomecov.bed > ${sampleid}_zero_coverage.bed
    bedtools maskfasta -fi ${vcf_applied_fasta} -bed ${sampleid}_zero_coverage.bed -fo ${sampleid}_bcftools_masked_consensus.fasta
    """
}
/*
process BBDUK { 
  tag "${sampleid}"
  label "setting_22"
  containerOptions "${bindOptions}"
  publishDir "${params.outdir}/${sampleid}/04_cleaned", mode: 'copy'

  input:
    tuple val(sampleid), path(fastq1), path(fastq2), path(ref)
  output:
    path("${sampleid}_rRNA_reads.log")
    path("${sampleid}_non_rRNA_fwd.fastq.gz")
    path("${sampleid}_non_rRNA_rev.fastq.gz")
    path("${sampleid}_rRNA_reads.log"), emit: bbduk_stats
    tuple val(sampleid), path("${sampleid}_non_rRNA_fwd.fastq.gz"), path("${sampleid}_non_rRNA_rev.fastq.gz"), emit: bbduk_filtered_fq

  script:
  """
  bbduk.sh -Xmx10g in=${fastq1} \
                   in2=${fastq2} \
                   out=${sampleid}_non_rRNA_fwd.fastq.gz \
                   out2=${sampleid}_non_rRNA_rev.fastq.gz \
                   outm=${sampleid}_rRNA_fwd.fastq.gz \
                   outm2=${sampleid}_rRNA_rev.fastq.gz \
                   k=31 ref=${ref} \
                   2>${sampleid}_rRNA_reads.log
  """
}

process FILTER_CONTROL { 
  tag "${sampleid}"
  label "setting_22"
  containerOptions "${bindOptions}"
  publishDir "${params.outdir}/${sampleid}/04_cleaned", mode: 'copy'

  input:
    tuple val(sampleid), path(fastq1), path(fastq2)
  output:
    path("${sampleid}_bbsplit_stats.txt")
    path("${sampleid}_cleaned_fwd.fastq.gz")
    path("${sampleid}_cleaned_rev.fastq.gz")
    path("${sampleid}_phyX_removed.fastq.gz")
    tuple val(sampleid), path("${sampleid}_bbsplit_stats.txt"), emit: stats
    path("${sampleid}_bbsplit_stats.txt"), emit: stats2
    tuple val(sampleid), path("${sampleid}_cleaned_fwd.fastq.gz"), path("${sampleid}_cleaned_rev.fastq.gz"), emit: bbsplit_filtered_fq

  script:
  """
  bbsplit.sh -Xmx10g ref=${params.phix} \
             in=${fastq1} \
             in2=${fastq2} \
             out_phiX174=${sampleid}_phyX_removed.fastq.gz \
             outu=${sampleid}_cleaned_fwd.fastq.gz \
             outu2=${sampleid}_cleaned_rev.fastq.gz \
             statsfile=${sampleid}_bbsplit_stats.txt
  """
}
*/
/*

*/
//Explore downtrack downloading krona taxonomy to see if it improves the visualisation
process KRAKEN2_TO_KRONA {
    tag "${sampleid}"
    label 'setting_3'
    publishDir "${params.outdir}/${sampleid}/05_read_classification", mode: 'copy'

    input:
    tuple val(sampleid), path(kraken_report)

    output:
    file("${sampleid}_kraken_krona.html")
    tuple val(sampleid), path("${sampleid}_kraken_krona.html")

    script:
    """
    ktImportText \\
        -o ${sampleid}_kraken_krona.html \\
        ${kraken_report}
    """
}


//the logic of the original est_abundance.py had to be modified as it was not working as intended for viral species 
// only defined at S1 (strain) level but not S level, these would just not appear in the bracken report. 
// Hence the updated_est_abundance.py script is used here instead.
//It will rescue the abundance estimates for such species by summing up all S1 level abundances to S level.
process BRACKEN {
  tag "${sampleid}"
  label 'setting_2'
  publishDir "$params.outdir/$sampleid/05_read_classification",  mode: 'copy'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(kraken_report)

  output:
    file("${sampleid}_bracken_report*.txt")
    tuple val(sampleid), path("${sampleid}_bracken_report.txt"), emit: bracken_results2
    tuple val(sampleid), path("${sampleid}_bracken_report_viral.txt"), emit: bracken_results

  script:
  """
  c1grep() { grep "\$@" || test \$? = 1; }

  updated_est_abundance.py -i ${kraken_report} \
                  -k ${params.kraken2_db}/database50mers.kmer_distrib \
                  -t 1 \
                  -l S -o ${sampleid}_bracken_report.txt


  c1grep  "taxonomy_id\\|virus\\|viroid" ${sampleid}_bracken_report.txt > ${sampleid}_bracken_report_viral.txt
  awk -F'\\t'  '\$7>=0.0001'  ${sampleid}_bracken_report_viral.txt > ${sampleid}_bracken_report_viral_filtered.txt
  """
}

//Explore downtrack downloading krona taxonomy to see if it improves the visualisation
process KRONA {
  publishDir "${params.outdir}/${sampleid}/05_read_classification", mode: 'link'
  label 'setting_3'
  containerOptions "${bindOptions}"
  tag "${sampleid}"

  input:
    tuple val(sampleid), path(krona_input)

  output:
    file "${sampleid}_kaiju_krona.html"

  script:
    """
    ktImportText -o ${sampleid}_kaiju_krona.html ${krona_input}
    """
}

process RETRIEVE_VIRAL_READS_KRAKEN2 {
  tag "${sampleid}"
  label "setting_10"
  containerOptions "${bindOptions}"
  publishDir "${params.outdir}/${sampleid}/05_read_classification", mode: 'copy'

  input:
    tuple val(sampleid), path(kraken_report), path(kraken_output), path(fastq1), path(fastq2), path(unc_fastq1), path(unc_fastq2)
  output:
    tuple val(sampleid), path("${sampleid}_cand_path_R1.fastq.gz"), path("${sampleid}_cand_path_R2.fastq.gz"), emit: fastq

  script:
  """
  extract_kraken_reads.py -k ${kraken_output} -r ${kraken_report} \
                          -t 10239 --include-children \
                          -s1 ${fastq1} -s2 ${fastq2} \
                          --fastq-output \
                          -o ${sampleid}_extracted_reads1.fastq -o2 ${sampleid}_extracted_reads2.fastq
  gzip ${sampleid}_extracted_reads1.fastq
  gzip ${sampleid}_extracted_reads2.fastq
  cat ${unc_fastq1} ${sampleid}_extracted_reads1.fastq.gz > ${sampleid}_cand_path_R1.fastq.gz
  cat ${unc_fastq2} ${sampleid}_extracted_reads2.fastq.gz >  ${sampleid}_cand_path_R2.fastq.gz
  """
}
/*
process DIAMOND  {
  tag "${sampleid}"
  label "setting_27"
  containerOptions "${bindOptions}"
  publishDir "${params.outdir}/${sampleid}/06_annotation", mode: 'copy'

  input:
    tuple val(sampleid), path(assembly)
  output:
    tuple val(sampleid), path("${sampleid}_diamond_matches.tsv"), emit: diamond_results

  script:
  """
  diamond blastx --query ${assembly} \
                 --db ${params.diamond_db} \
                 --out ${sampleid}_diamond_matches.tsv \
                 --outfmt 6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen evalue bitscore stitle staxids qseq sseq \
                 --evalue 1e-3 \
                 --max-target-seqs 5 \
                 --threads ${task.cpus}
  """
}
*/
process GENOMAD {
    tag "${sampleid}"
    label "setting_27"
    containerOptions "${bindOptions}"
    publishDir "${params.outdir}/${sampleid}/07_annotation/genomad", mode: 'copy'

    input:
      tuple val(sampleid), path(fasta)
    output:
      file "*_summary/*_virus.fna"
      file "*_summary/*_virus_summary.tsv"
      file "*_summary/*_virus_genes.tsv"
      file "*_summary/*_virus_proteins.faa"
      file "*_find_proviruses/*_provirus.tsv"
      file "*_find_proviruses/*_provirus_taxonomy.tsv"
      file "*_find_proviruses/*_provirus.fna"
      file "*_find_proviruses/*_provirus_genes.tsv"
      file "*_find_proviruses/*_provirus_proteins.faa"
      file "*_aggregated_classification/*aggregated_classification.tsv"
      file "*_annotate/*_taxonomy.tsv"
      file "${sampleid}_genomad.log"
      tuple val(sampleid), path("${sampleid}_scaffolds_virus_summary.tsv"), emit: virus_preds
      

    script:
    """
      genomad \\
        end-to-end \\
        ${fasta} \\
        ./ \\
        ${params.genomad_db} \\
        --threads ${task.cpus} \\
        --min-score 0.7 \\
        --splits 1 \\
        > ${sampleid}_genomad.log 2>&1
        cp ${sampleid}_scaffolds_summary/${sampleid}_scaffolds_virus_summary.tsv .
    """
}
//Derive ORFs from contig sequences using orfipy
//https://github.com/urmi-21/orfipy?tab=readme-ov-file
//Other options to consider are prodigal, OrfM and getorf

process ORFIPY {
  tag "${sampleid}"
    label "setting_21"
    containerOptions "${bindOptions}"
    publishDir "${params.outdir}/${sampleid}/07_annotation", mode: 'copy'

    input:
      tuple val(sampleid), path(fasta)
    output:
      file "${sampleid}_orfs.fasta"
      tuple val(sampleid), path("${sampleid}_orfs.fasta"), emit: orf_fasta

    script:
    """
      orfipy ${fasta} \\
        --outdir . \\
        --pep ${sampleid}_orfs.fasta \\
        --min 300 \\
        --procs ${task.cpus}
    """
}

process HMMSCAN {
  tag "${sampleid}"
    label "setting_20"
    containerOptions "${bindOptions}"
    publishDir "${params.outdir}/${sampleid}/07_annotation", mode: 'copy'

    input:
      tuple val(sampleid), path(fasta)
    output:
      file "${sampleid}_orfs.fasta"
      file "${sampleid}_hmmscan*_output.txt"
      tuple val(sampleid), path("${sampleid}_hmmscan_per_target_output.txt"), emit: hmmscan_preds

    script:
    """
      hmmscan --cpu ${task.cpus} \\
              --domtblout ${sampleid}_hmmscan_per_domain_output.txt \\
              --tblout ${sampleid}_hmmscan_per_target_output.txt \\
              --pfamtblout ${sampleid}_hmmscan_succinct_output.txt \\
              ${params.hmmer_db} ${fasta} \\
              > ${sampleid}_hmmscan.log 2>&1
    """
}

process SUMMARISE_RESULTS {
  tag "${sampleid}"
  label "setting_2"
  publishDir "${params.outdir}/${sampleid}/10_results_summary", mode: 'copy', pattern: '{*summary_viral_results.tsv}'
  publishDir "${params.outdir}/${sampleid}/07_annotation", mode: 'copy', pattern: '{*hmm_domain_summary_counts.tsv}'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(kraken_results), path(kaiju_results), path(blast), path(hmmscan)

  output:
    path("${sampleid}_summary_viral_results.tsv")
    path("${sampleid}_hmm_domain_summary_counts.tsv")
    tuple val(sampleid), path("${sampleid}_hmm_domain_summary_counts.tsv"), emit: domain_count
    tuple val(sampleid), path("${sampleid}_summary_viral_results.tsv"), emit: summary_known_viruses

  script:
  """
  viral_results_summary.py --kaiju ${kaiju_results} --sample_name ${sampleid} --kraken ${kraken_results} --blast ${blast} --hmmscan ${hmmscan}
  """
}

process NOVELS {
  tag "${sampleid}"
  label "setting_2"
  publishDir "${params.outdir}/${sampleid}/10_results_summary", mode: 'copy'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(contigs), path(hmmscan), path(genomad), path(blast)

  output:
    path("${sampleid}_novel_virus_candidates.tsv")
    tuple val(sampleid), path("${sampleid}_novel_virus_candidates.tsv"), emit: novel_virus_candidates


  script:
  """
  novel_candidates.py --sample_name ${sampleid} --fasta ${contigs} --genomad ${genomad} --hmmscan ${hmmscan} --blast ${blast}
  """
}

process EXTRACT_CONTIGS {
  tag "${sampleid}"
  label "setting_2"

  input:
    tuple val(sampleid), path(contig_ids), path(contigs)
    output:
    tuple val(sampleid), path("${sampleid}_candidate_viral_contigs.fasta"), emit: fasta

  script:
    """
    if [[ ! -s ${contig_ids} ]]; then
      touch ${sampleid}_candidate_viral_contigs.fasta
    else
      seqtk subseq ${contigs} ${contig_ids} > ${sampleid}_candidate_viral_contigs.fasta
    fi
    """
}
include { CAT_FASTQ } from './modules/cat_fastq/main'
include { FASTQC as FASTQC_RAW  } from './modules/fastqc/main'
include { FASTQC as FASTQC_TRIM } from './modules/fastqc/main'
include { FASTP } from './modules/fastp/main'
include { BBMAP_BBDUK } from './modules/bbmap/bbduk/main'
include { BBMAP_BBSPLIT } from './modules/bbmap/bbsplit/main'
include { KRAKEN2_KRAKEN2 } from './modules/kraken2/main'
include { KAIJU_KAIJU } from './modules/kaiju/main'

workflow {
  TIMESTAMP_START ()
  ch_versions = Channel.empty()
  /*
  if (params.samplesheet) {
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map { row ->
        // Check required fields
        if (!row.sampleid )  {
          exit 1, "ERROR: samplesheet is missing required fields for sample_id."
        }
        else if (!row.fastq_1 || !row.fastq_2) {
          exit 1, "ERROR: samplesheet is missing required field for one or both of the fastq files."
        }
        // Return parsed row
        tuple((row.sampleid), file(row.fastq_1), file(row.fastq_2)) }
      .set{ ch_sample }
  } else { exit 1, "Input samplesheet file not specified!" }
*/
  if (params.input) {
    Channel
          .fromSamplesheet("input")
          .map {
              meta, fastq_1, fastq_2 ->
                  if (!fastq_2) {
                      return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                  } else {
                      return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                  }
          }
          .groupTuple()
          .map { sample_id, meta_list, fastq_lists ->
              [meta_list[0], fastq_lists.flatten()]
          }
          .branch { meta, fastqs ->
              single  : fastqs.size() == 1
                  return [ meta, fastqs ]
              multiple: fastqs.size() > 1
                  return [ meta, fastqs ]
          }
          .set { ch_fastq }
  }
  //ch_fastq.single.view()
  //ch_fastq.multiple.view()
  CAT_FASTQ (
      ch_fastq.multiple
  )
  .reads
  .mix(ch_fastq.single)
  .set { ch_cat_fastq }
  ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

  configyaml = Channel.fromPath(workflow.commandLine.split(" -params-file ")[1].split(" ")[0])

  //FASTP ( ch_sample )
  FASTP ( CAT_FASTQ.out.reads, params.save_trimmed_fail, params.save_merged )
  trim_json         = FASTP.out.json
  trim_html         = FASTP.out.html
  trim_log          = FASTP.out.log
  trim_reads_fail   = FASTP.out.reads_fail
  trim_reads_merged = FASTP.out.reads_merged
  ch_versions       = ch_versions.mix(FASTP.out.versions.first())
  
  FASTQC_RAW ( CAT_FASTQ.out.reads )
  fastqc_raw_html = FASTQC_RAW.out.html
  fastqc_raw_zip  = FASTQC_RAW.out.zip
  ch_versions     = ch_versions.mix(FASTQC_RAW.out.versions.first())
  
  //Incorporate a logic like in nf-core/rnaseq where minimum reads has to be achieved to proceed?
  FASTQC_TRIM ( FASTP.out.reads )
  fastqc_trim_html = FASTQC_TRIM.out.html
  fastqc_trim_zip  = FASTQC_TRIM.out.zip
  ch_versions      = ch_versions.mix(FASTQC_TRIM.out.versions.first())
  
  //Filtering with sortmerna takes much longer than bbduk so use bbduk for prototype
  ch_rrna = Channel.fromPath(params.rrna_ref)
  BBMAP_BBDUK ( FASTP.out.reads, ch_rrna)
  //remove phiX reads
  bbmaplit_primary_ref_ch = Channel.fromPath(params.phix)
  empty_refs_ch = Channel.value( tuple([], []) )
  empty_index_ch         = Channel.empty()

  BBMAP_BBSPLIT ( BBMAP_BBDUK.out.reads, [], bbmaplit_primary_ref_ch, empty_refs_ch, false ) 
  
  //provide option to filter host or filter a plant host by default?

  //read classification with Kraken
  trial_ch = BBMAP_BBSPLIT.out.all_fastq.map { meta, reads ->
    def sample_id = meta.id
    def read1 = reads[0]
    def read2 = reads[1]
    tuple(sample_id, read1, read2)
  }
  stats_ch = BBMAP_BBSPLIT.out.stats.map { meta, stats ->
    def sample_id = meta.id
    def stats1 = stats
    tuple(sample_id, stats1)
  }

  KRAKEN2_KRAKEN2(BBMAP_BBSPLIT.out.all_fastq, params.kraken2_db, params.kraken2_save_classified_reads, params.kraken2_save_unclassified_reads, params.kraken2_save_readclassifications)
  BRACKEN ( KRAKEN2_KRAKEN2.out.report )
  //KRAKEN2_TO_KRONA ( KRAKEN2_KRAKEN2.out.kraken2_results2 )

  //retrieve reads that were not classified and reads classified as viral by kraken2
  //merge
  kraken_ch = KRAKEN2_KRAKEN2.out.results.map { meta, report, output, raw_reads, unclassified ->
    def sample_id = meta.id
    def read1 = raw_reads[0]
    def read2 = raw_reads[1]
    def unclassified1 = unclassified[0]
    def unclassified2 = unclassified[1]
    tuple(sample_id, report, output, read1, read2, unclassified1, unclassified2)
  }
  RETRIEVE_VIRAL_READS_KRAKEN2 ( kraken_ch )

  //read classification with kaiju
  //KAIJU ( FILTER_CONTROL.out.bbsplit_filtered_fq )
  ch_kaijudb = Channel.fromPath(params.kaiju_db_path)
  //incorporate a separate module for kaiju2krona and kaiju2table
  KAIJU_KAIJU (BBMAP_BBSPLIT.out.all_fastq, ch_kaijudb)
  //KAIJU_KAIJU (BBMAP_BBSPLIT.out.all_fastq)
  KRONA ( KAIJU_KAIJU.out.krona_results )
  read_classification_ch = KAIJU_KAIJU.out.kaiju_results.join(BRACKEN.out.bracken_results2)
                                                    .join(stats_ch)
  //                                                .join(FILTER_CONTROL.out.stats)

  SUMMARISE_READ_CLASSIFICATION ( read_classification_ch )

  //perform de novo assembly with spades using rnaspades
  SPADES ( RETRIEVE_VIRAL_READS_KRAKEN2.out.fastq )
  //Filter contigs by length less than 150 bp with SEQTK
  SEQTK ( SPADES.out.assembly )
  //Predict ORFs on filtered contigs
  ORFIPY ( SEQTK.out.filt_fasta )
  HMMSCAN ( ORFIPY.out.orf_fasta )

  GENOMAD ( SPADES.out.assembly )
  //Enhancement: Option to perform a blastx alignment of contig ORFs?
  //DIAMOND  ( SEQTK.out.filt_fasta.splitFasta(by: 5000, file: true) )
  //DIAMOND.out.diamond_results
  //  .groupTuple()
  //  .set { ch_blastxresults }

  BLASTN( SEQTK.out.filt_fasta.splitFasta(by: 5000, file: true) )
  BLASTN.out.blast_results
    .groupTuple()
    .set { ch_blastresults } 
  EXTRACT_VIRAL_BLAST_HITS ( ch_blastresults )
  //Add contig sequence to blast results summary table
  FASTA2TABLE ( EXTRACT_VIRAL_BLAST_HITS.out.viral_blast_results.join(SEQTK.out.filt_fasta) )
  //Mapping back to contigs that had viral blast hits
  EXTRACT_CONTIGS ( FASTA2TABLE.out.contig_ids.join(SEQTK.out.filt_fasta) )
  //MAPPING_BACK_TO_CONTIGS ( EXTRACT_CONTIGS.out.fasta.join(FILTER_CONTROL.out.bbsplit_filtered_fq) )
  MAPPING_BACK_TO_CONTIGS ( EXTRACT_CONTIGS.out.fasta.join(trial_ch) )
  SAMTOOLS_CONTIGS ( MAPPING_BACK_TO_CONTIGS.out.contig_aligned_sam )
  PYFAIDX_CONTIGS ( EXTRACT_CONTIGS.out.fasta )
  MOSDEPTH_CONTIGS (SAMTOOLS_CONTIGS.out.sorted_bam.join(PYFAIDX_CONTIGS.out.bed))

  contig_cov_stats_summary_ch = MOSDEPTH_CONTIGS.out.mosdepth_results.join(FASTA2TABLE.out.blast_results2)
                                                      .join(stats_ch)
                                                      .join(SAMTOOLS_CONTIGS.out.coverage)
                                                      .join(SAMTOOLS_CONTIGS.out.mapping_quality)
  CONTIG_COVSTATS(contig_cov_stats_summary_ch)
  //Mapping back to reference sequences retrieved from blast hits
  EXTRACT_REF_FASTA ( FASTA2TABLE.out.ref_ids )
  CLUSTER ( EXTRACT_REF_FASTA.out.fasta_files )
  //mapping_ch = CLUSTER.out.clusters.join(FILTER_CONTROL.out.bbsplit_filtered_fq)
  mapping_ch = CLUSTER.out.clusters.join(trial_ch)
  MAPPING_BACK_TO_REF ( mapping_ch )
  SAMTOOLS2 ( MAPPING_BACK_TO_REF.out.aligned_sam )
  BCFTOOLS ( SAMTOOLS2.out.sorted_bam )
  BEDTOOLS ( BCFTOOLS.out.vcf_applied_fasta )
  PYFAIDX ( EXTRACT_REF_FASTA.out.fasta_files )
  MOSDEPTH (SAMTOOLS2.out.sorted_bam.join(PYFAIDX.out.bed))
  cov_stats_summary_ch = MOSDEPTH.out.mosdepth_results.join(FASTA2TABLE.out.blast_results)
                                                      .join(stats_ch)
                                                      .join(BEDTOOLS.out.bcftools_masked_consensus_fasta)
                                                      .join(SAMTOOLS2.out.coverage)
                                                      .join(SAMTOOLS2.out.mapping_quality)
                                                      .join(EXTRACT_REF_FASTA.out.fasta_files)
  COVSTATS(cov_stats_summary_ch)
  FASTA2TABLE2  ( COVSTATS.out.detections_summary3)
  
  //Derive QC report
  // Merge all the  files into one channel
  //ch_multiqc_files = FASTP.out.fastp_json




  ch_multiqc_files = FASTP.out.json.map { meta, json ->
                      json
                      }
                      .mix(BBMAP_BBSPLIT.out.stats2)
                      .mix(BBMAP_BBDUK.out.log2)
                      .collect()

  QCREPORT(ch_multiqc_files)
  SUMMARISE_RESULTS ( SUMMARISE_READ_CLASSIFICATION.out.kraken_summary.join(SUMMARISE_READ_CLASSIFICATION.out.kaiju_summary)
                                                                      .join(CONTIG_COVSTATS.out.detections_summary) 
                                                                      .join(HMMSCAN.out.hmmscan_preds)
                                                                      )
  NOVELS ( SEQTK.out.filt_fasta.join(SUMMARISE_RESULTS.out.domain_count)
                    .join(GENOMAD.out.virus_preds)
                    .join(EXTRACT_VIRAL_BLAST_HITS.out.viral_blast_results)
         )
  fastqc_raw_html_fixed = fastqc_raw_html.map { meta, html ->
    tuple(meta.id, html)
  }

  fastqc_trim_html_fixed = fastqc_trim_html.map { meta, html ->
    tuple(meta.id, html)
  }

  trim_html_fixed = trim_html.map { meta, html ->
      tuple(meta.id, html)
  }

  
  files_for_report_ind_samples_ch = fastqc_raw_html_fixed.join(fastqc_trim_html_fixed)
                                                    .join(trim_html_fixed)
                                                    .join(SPADES.out.assembly)
                                                    .join(SUMMARISE_RESULTS.out.summary_known_viruses)                                               
                                                    .join(SUMMARISE_READ_CLASSIFICATION.out.kaiju_summary)
                                                    .join(SUMMARISE_READ_CLASSIFICATION.out.kraken_summary)
                                                    .join(FASTA2TABLE.out.blast_results2)
                                                    .join(FASTA2TABLE2.out.detections_summary_final)
                                                    .join(SAMTOOLS2.out.sorted_bam)
                                                    .join(NOVELS.out.novel_virus_candidates)
  
  //files_for_report_ind_samples_ch.view()
     
  files_for_report_global_ch = TIMESTAMP_START.out.timestamp
            .concat(QCREPORT.out.qc_report_html)
            .concat(QCREPORT.out.qc_report_txt)
            .concat(configyaml)
            .concat(Channel.from(params.input).map { file(it) }).toList()
  HTML_REPORT(files_for_report_ind_samples_ch
            .combine(files_for_report_global_ch))
  //files_for_report_global_ch.view()
/*
    if (params.subsample) {
      SUBSAMPLE ( REFORMAT.out.reformatted_fq )
      final_fq = SUBSAMPLE.out.subsampled_fq
    }
    else {
      final_fq = REFORMAT.out.reformatted_fq
    }
      //Perform direct alignment to a reference
      else if ( params.analysis_mode == 'map2ref') {
        MINIMAP2_REF ( final_fq )
        SAMTOOLS ( MINIMAP2_REF.out.aligned_sample )
        MEDAKA ( SAMTOOLS.out.sorted_sample )
        FILTER_VCF ( MEDAKA.out.unfilt_vcf )
      }

      else {
        error("Analysis mode (clustering) not specified with e.g. '--analysis_mode clustering' or via a detectable config file.")
      }
    }
  }
*/
}
