#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { fromSamplesheet } from 'plugin/nf-validation'

def helpMessage () {
    log.info """
    Virus Integrated Evaluation Workflow
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

def isNonEmptyFile(file) {
    return file.exists() && file.size() > 0
}

def buildBindOptions() {
  def bindOptions = "" 
  if (workflow.containerEngine == "singularity") {
    [ 
      params.blastn_db_dir,
      params.taxdump,
      params.kaiju_db_path,
      params.kraken2_db,
      params.genomad_db,
      params.hmmer_db_dir,
      params.prot_db
    ]
    .findAll { it } 
    .collect { file(it) }
    .collect { it.isDirectory() ? it : it.parent }
    .unique()
    .each { bindOptions += "-B ${it}:${it} " } 
  } 
  return bindOptions 
}

process EXTRACT_VIRAL_BLAST_HITS {
    tag "${sampleid}"
    label 'setting_7'
    containerOptions "--bind ${file(params.taxdump)}"

    input:
    tuple val(sampleid), path(blast_results), path(assembly_headers), path(fasta)
    path(taxonkit_db)

    output:
    tuple val(sampleid), path("${sampleid}_megablast_top_viral_hits.txt"), emit: viral_blast_results
    tuple val(sampleid), file("${sampleid}_viral_candidate_contig_ids.txt"), emit: viral_candidate_contig_ids
    tuple val(sampleid), file("${sampleid}_viral_candidate_contig_ids.txt"), file("${sampleid}_other_contig_ids.txt"), emit: contig_ids

    script:
    """
    cat ${blast_results} > ${sampleid}_blastn.txt
    filter_blast.py --blastn_results ${sampleid}_blastn.txt --sample_name ${sampleid} --taxonkit_database_dir ${taxonkit_db} --filter ${params.filter_terms} --assembly_headers ${assembly_headers}
    cut -f2 ${sampleid}_megablast_top_viral_hits.txt | sed '1d' | sed 's/ //g' | sort | uniq > ${sampleid}_viral_candidate_contig_ids.txt
    grep '>' ${fasta} | sed 's/>//g' | sed 's/ /_/g' | sort | uniq > ${sampleid}_all_contigs.txt
    grep -F -x -v -f ${sampleid}_viral_candidate_contig_ids.txt \
      ${sampleid}_all_contigs.txt \
      > ${sampleid}_other_contig_ids.txt || touch ${sampleid}_other_contig_ids.txt
    """
}

process EXTRACT_VIRAL_BLAST_HITS_ROUND2 {
    tag "${sampleid}"
    label 'setting_7'
    publishDir { "${params.outdir}/${sampleid}/07_annotation" }, mode: 'copy'

    input:
    tuple val(sampleid), path(blast_results), path(assembly_headers)
    path(taxonkit_db)

    output:
    tuple val(sampleid), path("${sampleid}_megablast_top_viral_hits.txt"), emit: viral_blast_results
    path("${sampleid}_blastn.txt")

    script:
    """
    cat ${blast_results} > ${sampleid}_blastn.txt
    filter_blast.py --blastn_results ${sampleid}_blastn.txt --sample_name ${sampleid} --taxonkit_database_dir ${taxonkit_db} --filter ${params.filter_terms} --assembly_headers ${assembly_headers}
    """
}

process MAPPING_BACK_TO_REF {
    tag "${sampleid}"
    label 'setting_10'

    input:
    tuple val(sampleid), path(ref), path(fastq1), path(fastq2)

    output:
    tuple val(sampleid), path(ref), file("${sampleid}_ref_aln.sam"), emit: aligned_sam

    script:
    """
    #bowtie2-build $ref $ref
    #bowtie2 --threads ${task.cpus} --very-sensitive-local -k 100 -x $ref \
    #        -1 $fastq1 -2 $fastq2 -S ${sampleid}_ref_aln.sam 2>> ${sampleid}_mapping.log
    bwa index ${ref}
    bwa mem -t ${task.cpus} ${ref} $fastq1 $fastq2 > ${sampleid}_ref_aln.sam 2>> ${sampleid}_mapping.log
    """
}

process REALIGN {
    tag "${sampleid}"
    label 'setting_10'

    input:
    tuple val(sampleid), path(contigs), path(fastq1), path(fastq2)

    output:
    tuple val(sampleid), path(contigs), file("${sampleid}_contig_realn.sam"), emit: contig_aligned_sam

    script:
    """
    bwa index ${contigs}
    bwa mem -t ${task.cpus} ${contigs} $fastq1 $fastq2 > ${sampleid}_contig_realn.sam 2>> ${sampleid}_mapping.log
    """
}

process SAMTOOLS2 {
    publishDir { "${params.outdir}/${sampleid}/09_mapping_to_ref" }, mode: 'copy'
    tag "${sampleid}"
    label 'setting_10'

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
    publishDir { "${params.outdir}/${sampleid}/08_mapping_to_contigs" }, mode: 'copy'
    tag "${sampleid}"
    label 'setting_10'

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

include { BBMAP_BBDUK } from './modules/bbmap/bbduk/main'
include { BBMAP_BBSPLIT } from './modules/bbmap/bbsplit/main'
include { BCFTOOLS } from './modules/bcftools/main'
include { BEDTOOLS } from './modules/bedtools/main'
include { BLAST_BLASTN as MEGABLAST } from './modules/blast/blastn/main'
include { BLAST_BLASTN as MEGABLAST_ROUND2 } from './modules/blast/blastn/main'
include { BLAST_BLASTN_TO_REF as MEGABLAST_TO_REF } from './modules/blast/blastn_to_ref/main'
include { BWA as MAP_TO_CONTIGS } from './modules/bwa/main'
include { CAT_FASTQ } from './modules/cat_fastq/main'
include { CDHIT_CDHIT as CLUSTER } from './modules/cdhit/cdhit/main'
include { COVSTATS as CONTIG_COVSTATS} from './modules/covstats/main'
include { COUNT_FASTQ_READS } from './modules/count_fastq_reads/main'
include { COVSTATS as REF_COVSTATS} from './modules/covstats/main'
include { DIAMOND_BLASTX } from './modules/diamond/blastx/main'
include { ENTREZDIRECT_EFETCH as EXTRACT_REF_FASTA } from './modules/entrezdirect/efetch/main'
include { EXTRACT_BLAST_HITS } from './modules/extract_blast_hits/main'
include { FASTA2TABLE_CONTIGS } from './modules/fasta2table_contigs/main'
include { FASTA2TABLE_REF } from './modules/fasta2table_ref/main'
include { FASTP } from './modules/fastp/main'
include { FASTQC as FASTQC_RAW } from './modules/fastqc/main'
include { FASTQC as FASTQC_TRIM } from './modules/fastqc/main'
include { FQ_SUBSAMPLE } from './modules/fq/subsample/main'
include { GENOMAD_ENDTOEND } from './modules/genomad/endtoend/main'
include { HMMSCAN } from './modules/hmmscan/main'
include { HTML_REPORT } from './modules/html_report/main'
include { IDENTIFY_ERRORS } from './modules/identify_errors/main'
include { KAIJU_KAIJU } from './modules/kaiju/main'
include { KRAKEN2_ABUNDANCE_ESTIMATE } from './modules/kraken2_abundance_estimate/main'
include { KRAKEN2_KRAKEN2 } from './modules/kraken2/main'
include { KRAKEN2_TO_KRONA } from './modules/kraken2_to_krona/main'
include { KRONA_KTIMPORTTEXT } from './modules/krona/ktimporttext/main'
include { MOSDEPTH as MOSDEPTH_CONTIGS } from './modules/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_REF } from './modules/mosdepth/main'
include { ORFIPY } from './modules/orfipy/main'
include { PYFAIDX as PYFAIDX_CONTIGS } from './modules/pyfaidx/main'
include { PYFAIDX as PYFAIDX_REF } from './modules/pyfaidx/main'
include { QC_REPORT } from './modules/qc_report/main'
include { RETRIEVE_VIRAL_READS_KRAKEN2 } from './modules/retrieve_viral_reads_kraken2/main'
include { SAMTOOLS_MPILEUP } from './modules/samtools/mpileup/main'
include { SEQTK_SAMPLE } from './modules/seqtk/sample/main'
include { SEQTK_SEQ } from './modules/seqtk/seq/main'
include { SEQTK_SUBSEQ as EXTRACT_CONTIGS } from './modules/seqtk/subseq/main'
include { SPADES } from './modules/spades/main'
include { START_TIMESTAMP } from './modules/start_timestamp/main'
include { SUMMARISE_RESULTS} from './modules/summarise_results/main'
include { SUMMARISE_READ_CLASSIFICATION } from './modules/summarise_read_classification/main'
include { TRIM_ENDS } from './modules/trim_ends/main'

workflow {
  // Show help message
  if (params.help) {
      helpMessage()
      exit 0
  }

  if (params.blastn_db) {
      blastn_db_name = file(params.blastn_db).name 
      params.blastn_db_dir = file(params.blastn_db).parent
  }

  if (params.hmmer_db) {
      params.hmmer_db_dir = file(params.hmmer_db).parent
  }
  
  if (params.prot_db) {
      prot_db_dir = file(params.prot_db).parent
  }

  params.bindOptions = buildBindOptions()
  //println "Bind options: ${params.bindOptions}"

  START_TIMESTAMP ()
  ch_versions = Channel.empty()
  
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
  CAT_FASTQ (
      ch_fastq.multiple
  )
  .reads
  .mix(ch_fastq.single)
  .set { ch_cat_fastq }
  ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

  configyaml = Channel.fromPath(workflow.commandLine.split(" -params-file ")[1].split(" ")[0])
  //configyaml = Channel.fromPath(params.config_yaml)
  //Probably best place to perform subsampling
  //Subsampling to 60M reads is as slow using the nf-core subsample module or seqtk sample (40-50 minutes)
    if ( params.subsample_enabled ) {
        //Check size of fastq file first before subsampling!
        //FQ_SUBSAMPLE ( BBMAP_BBSPLIT.out.all_fastq )
        ch_with_counts = CAT_FASTQ.out.reads \
            | COUNT_FASTQ_READS

        SEQTK_SAMPLE(
            ch_with_counts,
            params.subsample_size
        )

        ch_versions = ch_versions.mix( SEQTK_SAMPLE.out.versions.first() )
        merged_fastq = SEQTK_SAMPLE.out.reads.ifEmpty {
            CAT_FASTQ.out.reads
        }
    } else {
        merged_fastq = CAT_FASTQ.out.reads
    }

  /*
  //revisit, this logic was bugging but it would be ncie to combine COUNT_FASTQ_READS and SEQTK_SAMPLE into one process 
  //that performs subsampling if the file is above a certain size threshold, otherwise just passes through the original fastq file. This would avoid the issue of the channel not pairing after the first sample when using SEQTK_SAMPLE on its own.
  if ( params.subsample_enabled ) {
      //Check size of fastq file first before subsampling!
      //FQ_SUBSAMPLE ( BBMAP_BBSPLIT.out.all_fastq )
      //ch_with_counts = CAT_FASTQ.out.reads \
      //    | COUNT_FASTQ_READS

      //SEQTK_SAMPLE(
      //    ch_with_counts,
      //    params.subsample_size
      //)

      SEQTK_SAMPLE(
         CAT_FASTQ.out.reads,
         params.subsample_size
      )

      ch_versions = ch_versions.mix( SEQTK_SAMPLE.out.versions.first() )
      merged_fastq = SEQTK_SAMPLE.out.reads.ifEmpty {
          CAT_FASTQ.out.reads
      }
  } else {
      merged_fastq = CAT_FASTQ.out.reads
  }
  */

  FASTP ( merged_fastq, params.save_trimmed_fail, params.save_merged )
  trim_html         = FASTP.out.html
  trim_reads_for_fastqc   = FASTP.out.reads
  trim_reads_for_bbduk   = FASTP.out.reads
  ch_versions       = ch_versions.mix(FASTP.out.versions.first())
  
  FASTQC_RAW ( merged_fastq )
  fastqc_raw_html = FASTQC_RAW.out.html
  fastqc_raw_zip  = FASTQC_RAW.out.zip
  ch_versions     = ch_versions.mix(FASTQC_RAW.out.versions.first())
  
  FASTQC_TRIM ( trim_reads_for_fastqc )
  fastqc_trim_html = FASTQC_TRIM.out.html
  ch_versions      = ch_versions.mix(FASTQC_TRIM.out.versions.first())
  
  //Filtering with sortmerna takes much longer than bbduk so used bbduk for prototype
  //Nextflow zips channels together by default:
  //Task receives one item from trim_reads_for_bbduk + one item from ch_rrna
  //When using a chanmnel, for ex. ch_rrna, it has only one item, it stops pairing after the first sample, only one task runs
  //This is why only the first sample is processed.
  //Provide the rrna ref as a file parameter instead
  BBMAP_BBDUK ( trim_reads_for_bbduk, file(params.rrna_ref))

  //remove phiX reads
  BBMAP_BBSPLIT (
      BBMAP_BBDUK.out.reads,
      [],
      file(params.phix),
      [ [], [] ],
      false
  )

  //Provide option to filter host or filter a plant host by default?
  //READ CLASSIFICATION WITH KRAKEN2
  trial_ch = BBMAP_BBSPLIT.out.all_fastq.map { meta, reads ->
      def sample_id = meta.id
      def read1 = reads[0]
      def read2 = reads[1]
      tuple(sample_id, read1, read2)
  }
  stats_ch = BBMAP_BBSPLIT.out.stats.map { meta, stats ->
      def sampleid = meta.id
      def stats1 = stats
      tuple(sampleid, stats1)
  }
  KRAKEN2_KRAKEN2(BBMAP_BBSPLIT.out.all_fastq, params.kraken2_db, params.kraken2_save_classified_reads, params.kraken2_save_unclassified_reads, params.kraken2_save_readclassifications)
  
  //The logic of the original Bracken est_abundance.py had to be modified as it was not working as intended for viral species 
  // only defined at S1 (strain) level but not S level, these would just not appear in the bracken report. 
  //This script will rescue the abundance estimates for such species by summing up all S1 level abundances to S level.
  //It will also retain family and genus level that are not assigned.
  KRAKEN2_ABUNDANCE_ESTIMATE ( KRAKEN2_KRAKEN2.out.report )
  
  //explore downtrack downloading krona taxonomy to see if it improves the visualisation
  KRAKEN2_TO_KRONA ( KRAKEN2_KRAKEN2.out.report)

  //Retrieve reads that were not classified and reads classified as viral by kraken2
  //Merge
  kraken_ch = KRAKEN2_KRAKEN2.out.results.map { meta, report, output, raw_reads, unclassified ->
      def read1 = raw_reads[0]
      def read2 = raw_reads[1]
      def unclassified1 = unclassified[0]
      def unclassified2 = unclassified[1]
      tuple(meta, report, output, read1, read2, unclassified1, unclassified2)
  }

  RETRIEVE_VIRAL_READS_KRAKEN2 ( kraken_ch )

  //READ CLASSIFICATION WITH KAIJU
  //Incorporate a separate module for kaiju2krona and kaiju2table?
  //Explore downtrack downloading krona taxonomy to see if it improves the visualisation?
  KAIJU_KAIJU ( BBMAP_BBSPLIT.out.all_fastq, params.kaiju_db_path )
  KRONA_KTIMPORTTEXT ( KAIJU_KAIJU.out.krona_results )
  
  read_classification_ch = KAIJU_KAIJU.out.kaiju_results.join(KRAKEN2_ABUNDANCE_ESTIMATE.out.kraken2_results)
                                                    .join(stats_ch)
  SUMMARISE_READ_CLASSIFICATION ( read_classification_ch, params.taxdump )

  //perform de novo assembly with spades using rnaspades
  SPADES ( RETRIEVE_VIRAL_READS_KRAKEN2.out.fastq )
  //Filter contigs by length less than 150 bp with SEQTK
  SEQTK_SEQ ( SPADES.out.assembly )
  MEGABLAST( SEQTK_SEQ.out.filt_fasta.splitFasta(by: 2500, file: true), params.blastn_db )
  MEGABLAST.out.blast_results
      .groupTuple()
      .set { ch_blastresults }
  extract_viral_blast_hits_ch = ch_blastresults.join(SEQTK_SEQ.out.filt_headers)
      .join(SEQTK_SEQ.out.filt_fasta)
  EXTRACT_VIRAL_BLAST_HITS ( extract_viral_blast_hits_ch, params.taxdump )
  //Add contig sequence to blast results summary table
  //Mapping back to contigs that had viral blast hits
  EXTRACT_CONTIGS ( EXTRACT_VIRAL_BLAST_HITS.out.contig_ids.join(SEQTK_SEQ.out.filt_fasta) )
  MAP_TO_CONTIGS ( EXTRACT_CONTIGS.out.viral_candidate_fasta.join(trial_ch) )
  SAMTOOLS_MPILEUP ( MAP_TO_CONTIGS.out.contig_aligned_sam )
  IDENTIFY_ERRORS ( SAMTOOLS_MPILEUP.out.pileup )
  TRIM_ENDS ( IDENTIFY_ERRORS.out.trimmed_coords )
  REALIGN ( TRIM_ENDS.out.trimmed_contigs.join(trial_ch) )
  SAMTOOLS_CONTIGS ( REALIGN.out.contig_aligned_sam )
  pyfaidx_contigs_input_ch = TRIM_ENDS.out.trimmed_contigs.map { sampleid, fasta ->
          tuple(sampleid, fasta, 'contigs')
  }
  pyfaidx_contigs = PYFAIDX_CONTIGS ( pyfaidx_contigs_input_ch )
  MOSDEPTH_CONTIGS (SAMTOOLS_CONTIGS.out.sorted_bam.join(pyfaidx_contigs.bed))
  MEGABLAST_ROUND2 ( TRIM_ENDS.out.trimmed_contigs, file(params.blastn_db) )
  extract_viral_blast_hits_round2_ch = MEGABLAST_ROUND2.out.blast_results.join(SEQTK_SEQ.out.filt_headers) 
  EXTRACT_VIRAL_BLAST_HITS_ROUND2 ( extract_viral_blast_hits_round2_ch, params.taxdump )
  fasta2table_contigs_input_ch = EXTRACT_VIRAL_BLAST_HITS_ROUND2.out.viral_blast_results
      .join(TRIM_ENDS.out.trimmed_contigs)
      .map { sampleid, tophits, fasta -> tuple(sampleid, tophits, fasta, 'contigs') }
  fasta2table_contigs = FASTA2TABLE_CONTIGS ( fasta2table_contigs_input_ch )
  contig_cov_stats_summary_ch = MOSDEPTH_CONTIGS.out.mosdepth_results.join(fasta2table_contigs.blast_results2)
      .join(stats_ch)
      .join(SAMTOOLS_CONTIGS.out.coverage)
      .join(SAMTOOLS_CONTIGS.out.mapping_quality)
      .map { sampleid, bed, blast_results, bbsplit_stats, coverage, mapping_q
      -> tuple(sampleid, 'contig', bed, blast_results, bbsplit_stats, coverage, mapping_q) }

  //Predict ORFs on filtered contigs
  //Derive ORFs from contig sequences using orfipy
  //https://github.com/urmi-21/orfipy?tab=readme-ov-file
  //Other options to consider are prodigal, OrfM and getorf 
  ORFIPY ( TRIM_ENDS.out.trimmed_contigs.join(EXTRACT_CONTIGS.out.other_fasta) )
  HMMSCAN ( ORFIPY.out.orf_fasta, params.hmmer_db )
  genomad_ch = TRIM_ENDS.out.trimmed_contigs.join(EXTRACT_CONTIGS.out.other_fasta)
  GENOMAD_ENDTOEND ( genomad_ch, params.genomad_db )

  //Enhancement: Option to perform a blastx alignment of contig ORFs?
  diamond_ch = TRIM_ENDS.out.trimmed_contigs.join(EXTRACT_CONTIGS.out.other_fasta) 
  DIAMOND_BLASTX ( diamond_ch, params.prot_db )
  CONTIG_COVSTATS( contig_cov_stats_summary_ch)
  //Mapping back to reference sequences retrieved from blast hits
  EXTRACT_REF_FASTA ( fasta2table_contigs.ref_ids )
  CLUSTER ( EXTRACT_REF_FASTA.out.fasta_files )
  MEGABLAST_TO_REF ( fasta2table_contigs.contig_fasta.join(CLUSTER.out.clusters) )
  mapping_ch = CLUSTER.out.clusters.join(trial_ch)
  MAPPING_BACK_TO_REF ( mapping_ch )
  SAMTOOLS2 ( MAPPING_BACK_TO_REF.out.aligned_sam )
  BCFTOOLS ( SAMTOOLS2.out.sorted_bam )
  BEDTOOLS ( BCFTOOLS.out.vcf_applied_fasta )
  pyfaidx_ref_input_ch = EXTRACT_REF_FASTA.out.fasta_files.map { sampleid, fasta ->
          tuple(sampleid, fasta, 'reference')
  }
  pyfaidx_ref = PYFAIDX_REF ( pyfaidx_ref_input_ch )
  MOSDEPTH_REF (SAMTOOLS2.out.sorted_bam.join(pyfaidx_ref.bed))
  ref_cov_stats_summary_ch =  MOSDEPTH_REF.out.mosdepth_results.join(fasta2table_contigs.blast_results)
      .join(stats_ch)
      .join(SAMTOOLS2.out.coverage)
      .join(SAMTOOLS2.out.mapping_quality)
      .map { sampleid, bed, blast_results, bbsplit_stats, coverage, mapping_q
      -> tuple(sampleid, 'reference', bed, blast_results, bbsplit_stats, coverage, mapping_q) }

  REF_COVSTATS(ref_cov_stats_summary_ch)
  fasta2table_ref_input_ch = REF_COVSTATS.out.detections_summary
      .join(BEDTOOLS.out.bcftools_masked_consensus_fasta)
      .map { sampleid, stats, fasta -> tuple(sampleid, stats, fasta, 'reference') }
  fasta2table_ref = FASTA2TABLE_REF ( fasta2table_ref_input_ch )
  
  //Derive QC report
  //Merge all the  files into one channel
  ch_multiqc_files = FASTP.out.json.map { meta, json ->
      json
      }
      .mix(BBMAP_BBSPLIT.out.stats2)
      .mix(BBMAP_BBDUK.out.log2)
      .collect()

  QC_REPORT(ch_multiqc_files)
  summarise_results_input_ch = SUMMARISE_READ_CLASSIFICATION.out.kraken_summary
      .join(SUMMARISE_READ_CLASSIFICATION.out.kaiju_summary)
      .join(CONTIG_COVSTATS.out.detections_summary)
      .join(HMMSCAN.out.hmmscan_preds)
      .join(fasta2table_ref.detections_summary_final)
      .join(SEQTK_SEQ.out.filt_fasta)
      .join(GENOMAD_ENDTOEND.out.virus_preds)
      .join(EXTRACT_VIRAL_BLAST_HITS.out.viral_blast_results)
      .join(DIAMOND_BLASTX.out.diamond_results)
      .map { sampleid, kraken_results, kaiju_results, blast, hmmscan, map2ref, contigs, genomad, blast_novel, diamond_results ->
        tuple(sampleid, kraken_results, kaiju_results, blast, hmmscan, map2ref, contigs, genomad, blast_novel, diamond_results, file(params.rvdb_taxonomy))
      }

  SUMMARISE_RESULTS(summarise_results_input_ch)
  
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
      .join(CONTIG_COVSTATS.out.detections_summary)
      .join(fasta2table_ref.detections_summary_final)
      .join(SAMTOOLS2.out.sorted_bam)
      .join(SUMMARISE_RESULTS.out.novel_virus_candidates)
      .join(SUMMARISE_RESULTS.out.novel_support_summary)
      .join(MEGABLAST_TO_REF.out.blast_results)
      .join((ORFIPY.out.orf_fasta))
      .join(HMMSCAN.out.hmmscan_domain_preds)
     
  files_for_report_global_ch = START_TIMESTAMP.out.timestamp
      .concat(QC_REPORT.out.qc_report_html)
      .concat(QC_REPORT.out.qc_report_txt)
      .concat(configyaml)
      .concat(Channel.from(params.input).map { file(it) }).toList()
  HTML_REPORT(files_for_report_ind_samples_ch
      .combine(files_for_report_global_ch))
}