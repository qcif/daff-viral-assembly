nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { fromSamplesheet } from 'plugin/nf-validation'
def isNonEmptyFile(file) {
    return file.exists() && file.size() > 0
}

include { BBMAP_BBDUK } from '../../modules/bbmap/bbduk/main'
include { BBMAP_BBSPLIT } from '../../modules/bbmap/bbsplit/main'
include { BCFTOOLS } from '../../modules/bcftools/main'
include { BEDTOOLS } from '../../modules/bedtools/main'
include { BLAST_BLASTN as MEGABLAST } from '../../modules/blast/blastn/main'
include { BLAST_BLASTN as MEGABLAST_ROUND2 } from '../../modules/blast/blastn/main'
include { BLAST_BLASTN_TO_REF as MEGABLAST_TO_REF } from '../../modules/blast/blastn_to_ref/main'
include { BWA_MEM as MAP_TO_CONTIGS } from '../../modules/bwa_mem/main'
include { BWA_MEM as MAPPING_BACK_TO_REF } from '../../modules/bwa_mem/main'
include { BWA_MEM as REALIGN } from '../../modules/bwa_mem/main'
include { CAT_FASTQ } from '../../modules/cat_fastq/main'
include { CDHIT_CDHIT as CLUSTER } from '../../modules/cdhit/cdhit/main'
include { COVSTATS as CONTIG_COVSTATS} from '../../modules/covstats/main'
include { COUNT_FASTQ_READS } from '../../modules/count_fastq_reads/main'
include { COVSTATS as REF_COVSTATS} from '../../modules/covstats/main'
include { DIAMOND_BLASTX } from '../../modules/diamond/blastx/main'
include { ENTREZDIRECT_EFETCH as EXTRACT_REF_FASTA } from '../../modules/entrezdirect/efetch/main'
include { EXTRACT_BLAST_HITS } from '../../modules/extract_blast_hits/main'
include { EXTRACT_RAW_VIRAL_BLAST_HITS } from '../../modules/extract_raw_viral_blast_hits/main'
include { EXTRACT_FINAL_VIRAL_BLAST_HITS } from '../../modules/extract_final_viral_blast_hits/main'
include { FASTA2TABLE_CONTIGS } from '../../modules/fasta2table_contigs/main'
include { FASTA2TABLE_REF } from '../../modules/fasta2table_ref/main'
include { FASTP } from '../../modules/fastp/main'
include { FASTQC as FASTQC_RAW } from '../../modules/fastqc/main'
include { FASTQC as FASTQC_TRIM } from '../../modules/fastqc/main'
include { FQ_SUBSAMPLE } from '../../modules/fq/subsample/main'
include { GENOMAD_ENDTOEND } from '../../modules/genomad/endtoend/main'
include { GENOMAD_DOWNLOAD_DB } from '../../modules/genomad/download_db/main'
include { HMMSCAN } from '../../modules/hmmscan/main'
include { HTML_REPORT } from '../../modules/html_report/main'
include { IDENTIFY_ERRORS } from '../../modules/identify_errors/main'
include { KAIJU_KAIJU } from '../../modules/kaiju/main'
include { KRAKEN2_ABUNDANCE_ESTIMATE } from '../../modules/kraken2_abundance_estimate/main'
include { KRAKEN2_KRAKEN2 } from '../../modules/kraken2/main'
include { KRAKEN2_TO_KRONA } from '../../modules/kraken2_to_krona/main'
include { KRONA_KTIMPORTTEXT } from '../../modules/krona/ktimporttext/main'
include { MOSDEPTH as MOSDEPTH_CONTIGS } from '../../modules/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_REF } from '../../modules/mosdepth/main'
include { ORFIPY } from '../../modules/orfipy/main'
include { PYFAIDX as PYFAIDX_CONTIGS } from '../../modules/pyfaidx/main'
include { PYFAIDX as PYFAIDX_REF } from '../../modules/pyfaidx/main'
include { QC_REPORT } from '../../modules/qc_report/main'
include { RETRIEVE_VIRAL_READS_KRAKEN2 } from '../../modules/retrieve_viral_reads_kraken2/main'
include { SAMTOOLS_CONTIGS } from '../../modules/samtools/contigs/main'
include { SAMTOOLS_MPILEUP } from '../../modules/samtools/mpileup/main'
include { SAMTOOLS_REF } from '../../modules/samtools/ref/main'
include { SEQTK_SAMPLE } from '../../modules/seqtk/sample/main'
include { SEQTK_SEQ } from '../../modules/seqtk/seq/main'
include { SEQTK_SUBSEQ as EXTRACT_CONTIGS } from '../../modules/seqtk/subseq/main'
include { SPADES } from '../../modules/spades/main'
include { START_TIMESTAMP } from '../../modules/start_timestamp/main'
include { SUMMARISE_RESULTS} from '../../modules/summarise_results/main'
include { SUMMARISE_READ_CLASSIFICATION } from '../../modules/summarise_read_classification/main'
include { TRIM_ENDS } from '../../modules/trim_ends/main'

workflow VIEW {
    // Show help message
    
    if ( !params.taxdump ) {
        error "Required parameter 'taxdump' is missing. Please set it in your -params-file."
    }
    else {
        params.taxdump_dir = file(params.taxdump).parent
    }

    if ( !params.kaiju_db ) {
        error "Required parameter 'kaiju_db' is missing. Please set it in your -params-file."
    }
    else {
        params.kaiju_db_dir = file(params.kaiju_db).parent
    }
    if ( !params.genomad_db) {
        if (workflow.profile.tokenize(',').contains('test')) {
            db_results = GENOMAD_DOWNLOAD_DB()
        }
        else {
            error "Required parameter 'genomad_db' is missing. Please set it in your -params-file."
        }
        ch_genomad_db = db_results.db
    }
    else {
        ch_genomad_db = Channel.fromPath(params.genomad_db)
    }
    
    def otherRequiredParams = [
        'blastn_db',
        'hmmer_db',
        'taxdump',
        'prot_db',
        'rvdb_taxonomy',
        'rrna_ref',
        'kraken2_db',
    ]

    otherRequiredParams.each { p ->
        if (!params[p]) {
            error "Required parameter '${p}' is missing. Please set it in your -params-file."
        }
    }
    

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

    //This is a bit hacky but it allows us to capture the params file used in the run and pass it to the report module without having to specify it as an output in every process. 
    //It also allows us to use a default params file for testing when the user does not specify one.
    def yamlFile

    if (workflow.commandLine.contains('-params-file')) {
        yamlFile = workflow.commandLine.split(" -params-file ")[1].split(" ")[0]
    }
    else if (workflow.profile.tokenize(',').contains('test')) {
        yamlFile = "${projectDir}/params/user_params_test.yml"
    }
    
    configyaml = Channel.fromPath(yamlFile)
    
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
    //revisit, this logic was bugging but it would be nice to combine COUNT_FASTQ_READS and SEQTK_SAMPLE into one process 
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
    ch_cleaned_fq = BBMAP_BBSPLIT.out.all_fastq.map { meta, reads ->
        def sample_id = meta.id
        def read1 = reads[0]
        def read2 = reads[1]
        tuple(sample_id, read1, read2)
    }
    ch_stats = BBMAP_BBSPLIT.out.stats.map { meta, stats ->
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
    ch_kraken = KRAKEN2_KRAKEN2.out.results.map { meta, report, output, raw_reads, unclassified ->
        def read1 = raw_reads[0]
        def read2 = raw_reads[1]
        def unclassified1 = unclassified[0]
        def unclassified2 = unclassified[1]
        tuple(meta, report, output, read1, read2, unclassified1, unclassified2)
    }

    RETRIEVE_VIRAL_READS_KRAKEN2 ( ch_kraken )

    //READ CLASSIFICATION WITH KAIJU
    //Incorporate a separate module for kaiju2krona and kaiju2table?
    //Explore downtrack downloading krona taxonomy to see if it improves the visualisation?
    KAIJU_KAIJU ( BBMAP_BBSPLIT.out.all_fastq, params.kaiju_db_dir )
    KRONA_KTIMPORTTEXT ( KAIJU_KAIJU.out.krona_results )
    
    ch_read_classification = KAIJU_KAIJU.out.kaiju_results.join(KRAKEN2_ABUNDANCE_ESTIMATE.out.kraken2_results)
                                                        .join(ch_stats)
    SUMMARISE_READ_CLASSIFICATION ( ch_read_classification, params.taxdump )

    //perform de novo assembly with spades using rnaspades
    SPADES ( RETRIEVE_VIRAL_READS_KRAKEN2.out.fastq )
    //Filter contigs by length less than 150 bp with SEQTK
    SEQTK_SEQ ( SPADES.out.assembly )
    ch_blast_db = Channel.value(
        tuple(
            file(params.blastn_db).parent,
            file(params.blastn_db).name
        )
    )

    MEGABLAST(
        SEQTK_SEQ.out.filt_fasta.splitFasta(by: 2500, file: true),
        ch_blast_db
    )
    MEGABLAST.out.blast_results
        .groupTuple()   
        .set { ch_blastresults }
    ch_extract_raw_viral_blast_hits = ch_blastresults
        .join(SEQTK_SEQ.out.filt_headers)
    EXTRACT_RAW_VIRAL_BLAST_HITS(ch_extract_raw_viral_blast_hits, params.taxdump)
    //Add contig sequence to blast results summary table
    //Mapping back to contigs that had viral blast hits
    EXTRACT_CONTIGS ( EXTRACT_RAW_VIRAL_BLAST_HITS.out.viral_blast_results.join(SEQTK_SEQ.out.filt_fasta) )
    MAP_TO_CONTIGS(
        EXTRACT_CONTIGS.out.viral_candidate_fasta
            .join(ch_cleaned_fq)
            .map { sampleid, contigs, fastq1, fastq2 ->
                tuple(sampleid, "contig_aln", contigs, fastq1, fastq2)
            }
    )

    SAMTOOLS_MPILEUP ( MAP_TO_CONTIGS.out.aligned_sam )
    IDENTIFY_ERRORS ( SAMTOOLS_MPILEUP.out.pileup )
    TRIM_ENDS ( IDENTIFY_ERRORS.out.trimmed_coords )
    REALIGN ( 
        TRIM_ENDS.out.trimmed_contigs.join(ch_cleaned_fq)
        .map { sampleid, contigs, fastq1, fastq2 ->
            tuple(sampleid, "contig_realn", contigs, fastq1, fastq2)
        }
    )
    SAMTOOLS_CONTIGS ( REALIGN.out.aligned_sam )
    ch_pyfaidx_contigs_input = TRIM_ENDS.out.trimmed_contigs.map { sampleid, fasta ->
            tuple(sampleid, fasta, 'contigs')
    }
    pyfaidx_contigs = PYFAIDX_CONTIGS ( ch_pyfaidx_contigs_input )
    MOSDEPTH_CONTIGS (SAMTOOLS_CONTIGS.out.sorted_bam.join(pyfaidx_contigs.bed))
    MEGABLAST_ROUND2 ( TRIM_ENDS.out.trimmed_contigs, ch_blast_db )
    ch_extract_final_viral_blast_hits = MEGABLAST_ROUND2.out.blast_results.join(SEQTK_SEQ.out.filt_headers) 
    EXTRACT_FINAL_VIRAL_BLAST_HITS ( ch_extract_final_viral_blast_hits, params.taxdump )
    ch_fasta2table_contigs_input = EXTRACT_FINAL_VIRAL_BLAST_HITS.out.viral_blast_results
        .join(TRIM_ENDS.out.trimmed_contigs)
        .map { sampleid, tophits, fasta -> tuple(sampleid, tophits, fasta, 'contigs') }
    ch_fasta2table_contigs = FASTA2TABLE_CONTIGS ( ch_fasta2table_contigs_input )
    ch_contig_cov_stats_summary = MOSDEPTH_CONTIGS.out.mosdepth_results.join(ch_fasta2table_contigs.blast_results2)
        .join(ch_stats)
        .join(SAMTOOLS_CONTIGS.out.coverage)
        .join(SAMTOOLS_CONTIGS.out.mapping_quality)
        .map { sampleid, bed, blast_results, bbsplit_stats, coverage, mapping_q
        -> tuple(sampleid, 'contig', bed, blast_results, bbsplit_stats, coverage, mapping_q) }

    //Predict ORFs on filtered contigs
    //Derive ORFs from contig sequences using orfipy
    //https://github.com/urmi-21/orfipy?tab=readme-ov-file
    //Other options to consider are prodigal, OrfM and getorf 
    ORFIPY ( TRIM_ENDS.out.trimmed_contigs.join(EXTRACT_CONTIGS.out.other_fasta) )
    ch_hmmer_db = Channel.value(
        tuple(
            file(params.hmmer_db).parent,
            file(params.hmmer_db).name
        )
    )

    HMMSCAN ( ORFIPY.out.orf_fasta, ch_hmmer_db )
    ch_genomad = TRIM_ENDS.out.trimmed_contigs.join(EXTRACT_CONTIGS.out.other_fasta)
    //GENOMAD_ENDTOEND ( genomad_ch, params.genomad_db )
    GENOMAD_ENDTOEND ( ch_genomad, ch_genomad_db )

    //Enhancement: Option to perform a blastx alignment of contig ORFs?
    ch_diamond = TRIM_ENDS.out.trimmed_contigs.join(EXTRACT_CONTIGS.out.other_fasta) 
    DIAMOND_BLASTX ( ch_diamond, params.prot_db )
    CONTIG_COVSTATS( ch_contig_cov_stats_summary)
    //Mapping back to reference sequences retrieved from blast hits
    EXTRACT_REF_FASTA ( ch_fasta2table_contigs.ref_ids )
    CLUSTER ( EXTRACT_REF_FASTA.out.fasta_files )
    MEGABLAST_TO_REF ( ch_fasta2table_contigs.contig_fasta.join(CLUSTER.out.clusters) )
    ch_mapping = CLUSTER.out.clusters.join(ch_cleaned_fq) 
        .map { sampleid, cluster_fasta, fastq1, fastq2 ->
            tuple(sampleid, "ref_aln", cluster_fasta, fastq1, fastq2)
        }

    MAPPING_BACK_TO_REF ( ch_mapping )
    SAMTOOLS_REF ( MAPPING_BACK_TO_REF.out.aligned_sam )
    BCFTOOLS ( SAMTOOLS_REF.out.sorted_bam )
    BEDTOOLS ( BCFTOOLS.out.vcf_applied_fasta )
    ch_pyfaidx_ref_input = EXTRACT_REF_FASTA.out.fasta_files.map { sampleid, fasta ->
            tuple(sampleid, fasta, 'reference')
    }
    ch_pyfaidx_ref = PYFAIDX_REF ( ch_pyfaidx_ref_input )
    MOSDEPTH_REF (SAMTOOLS_REF.out.sorted_bam.join(ch_pyfaidx_ref.bed))
    ch_ref_cov_stats_summary =  MOSDEPTH_REF.out.mosdepth_results.join(ch_fasta2table_contigs.blast_results)
        .join(ch_stats)
        .join(SAMTOOLS_REF.out.coverage)
        .join(SAMTOOLS_REF.out.mapping_quality)
        .map { sampleid, bed, blast_results, bbsplit_stats, coverage, mapping_q
        -> tuple(sampleid, 'reference', bed, blast_results, bbsplit_stats, coverage, mapping_q) }

    REF_COVSTATS(ch_ref_cov_stats_summary)
    ch_fasta2table_ref_input = REF_COVSTATS.out.detections_summary
        .join(BEDTOOLS.out.bcftools_masked_consensus_fasta)
        .map { sampleid, stats, fasta -> tuple(sampleid, stats, fasta, 'reference') }
    ch_fasta2table_ref = FASTA2TABLE_REF ( ch_fasta2table_ref_input )
    
    //Derive QC report
    //Merge all the  files into one channel
    ch_multiqc_files = FASTP.out.json.map { meta, json ->
        json
        }
        .mix(BBMAP_BBSPLIT.out.stats2)
        .mix(BBMAP_BBDUK.out.log2)
        .collect()

    QC_REPORT(ch_multiqc_files)
    ch_summarise_results_input = SUMMARISE_READ_CLASSIFICATION.out.kraken_summary
        .join(SUMMARISE_READ_CLASSIFICATION.out.kaiju_summary)
        .join(CONTIG_COVSTATS.out.detections_summary)
        .join(HMMSCAN.out.hmmscan_preds)
        .join(ch_fasta2table_ref.detections_summary_final)
        .join(SEQTK_SEQ.out.filt_fasta)
        .join(GENOMAD_ENDTOEND.out.virus_preds)
        .join(EXTRACT_FINAL_VIRAL_BLAST_HITS.out.viral_blast_results)
        .join(DIAMOND_BLASTX.out.diamond_results)
        .map { sampleid, kraken_results, kaiju_results, blast, hmmscan, map2ref, contigs, genomad, blast_novel, diamond_results ->
            tuple(sampleid, kraken_results, kaiju_results, blast, hmmscan, map2ref, contigs, genomad, blast_novel, diamond_results, file(params.rvdb_taxonomy))
        }

    SUMMARISE_RESULTS(ch_summarise_results_input)
    
    fastqc_raw_html_fixed = fastqc_raw_html.map { meta, html ->
        tuple(meta.id, html)
    }

    fastqc_trim_html_fixed = fastqc_trim_html.map { meta, html ->
        tuple(meta.id, html)
    }

    trim_html_fixed = trim_html.map { meta, html ->
        tuple(meta.id, html)
    }

    ch_files_for_report_ind_samples = fastqc_raw_html_fixed.join(fastqc_trim_html_fixed)
        .join(trim_html_fixed)
        .join(SPADES.out.assembly)
        .join(SUMMARISE_RESULTS.out.summary_known_viruses)                                               
        .join(SUMMARISE_READ_CLASSIFICATION.out.kaiju_summary)
        .join(SUMMARISE_READ_CLASSIFICATION.out.kraken_summary)
        .join(CONTIG_COVSTATS.out.detections_summary)
        .join(ch_fasta2table_ref.detections_summary_final)
        .join(SAMTOOLS_REF.out.sorted_bam)
        .join(SUMMARISE_RESULTS.out.novel_support_summary)
        .join(MEGABLAST_TO_REF.out.blast_results)
        .join((ORFIPY.out.orf_fasta))
        .join(HMMSCAN.out.hmmscan_domain_preds)
        .join(SUMMARISE_RESULTS.out.diamond_summary)
        .join(SUMMARISE_RESULTS.out.novel_contig_summary)
        
    ch_files_for_report_global = START_TIMESTAMP.out.timestamp
        .concat(QC_REPORT.out.qc_report_html)
        .concat(QC_REPORT.out.qc_report_txt)
        .concat(configyaml)
        .concat(Channel.from(params.input).map { file(it) }).toList()
    HTML_REPORT(ch_files_for_report_ind_samples
        .combine(ch_files_for_report_global))
}