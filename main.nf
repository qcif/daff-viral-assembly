#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

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

      #### Pre-processing and QC options ####
      --merge                         Merge fastq files with the same sample name
                                      Default: true
      --qc_only                       Only perform preliminary QC step using Nanoplot
                                      Default: false
      --preprocessing_only            Only perform preprocessing steps specied
                                      Default: false
      --porechop_options              Porechop_ABI options
                                      Default: ''
      --porechop_custom_primers       Limit porechop search to custom adapters specified under porechop_custom_primers_path
                                      Default: ''
      --porechop_custom_primers_path  Path to custom adpaters for porechop
                                      Default: ''
      --qual_filt                     Run quality filtering step using chopper
                                      [False]
      --chopper_options               Chopper options
                                      Default: ''

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
if (params.blastn_COI != null) {
    blastn_COI_name = file(params.blastn_COI).name
    blastn_COI_dir = file(params.blastn_COI).parent
}

if (params.sortmerna_ref != null) {
    sortmerna_ref_name = file(params.sortmerna_ref).name
    sortmerna_ref_dir = file(params.sortmerna_ref).parent
}

if (params.kraken2_db != null) {
    krkdb_dir = file(params.kraken2_db).parent
}
if (params.kaiju_nodes != null & params.kaiju_dbname != null & params.kaiju_names != null) {
    kaiju_dbs_dir = file(params.kaiju_nodes).parent
}
if (params.diamond_db != null) {
    diamond_db_dir = file(params.diamond_db).parent
}

if (params.genomad_db != null) {
    genomad_db_dir = file(params.genomad_db).parent
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

if (params.porechop_custom_primers == true) {
    porechop_custom_primers_dir = file(params.porechop_custom_primers_path).parent
}

def isNonEmptyFile(file) {
    return file.exists() && file.size() > 0
}

switch (workflow.containerEngine) {
  case "singularity":
    bindbuild = "";
    if (params.blastn_db != null) {
      bindbuild = (bindbuild + "-B ${blastn_db_dir} ")
    }
    if (params.blastn_COI != null) {
      bindbuild = (bindbuild + "-B ${blastn_COI_dir} ")
    }
    if (params.taxdump != null) {
      bindbuild = (bindbuild + "-B ${params.taxdump} ")
    }
    if (params.diamond_db != null) {
      bindbuild = (bindbuild + "-B ${diamond_db_dir} ")
    }

    if (params.sortmerna_ref != null) {
      bindbuild = (bindbuild + "-B ${sortmerna_ref_dir} ")
    }
    if (params.kraken2_db != null) {
      bindbuild = (bindbuild + "-B ${krkdb_dir} ")
    }
    if (params.kaiju_nodes != null & params.kaiju_dbname != null & params.kaiju_names != null) {
      bindbuild = (bindbuild + "-B ${kaiju_dbs_dir} ")
    }
    if (params.genomad_db != null ) {
      bindbuild = (bindbuild + "-B ${genomad_db_dir} ")
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

process FASTP {
  tag "$sampleid"
  publishDir "${params.outdir}/${sampleid}/02_qtrimmed", mode: 'copy', pattern: '{*fastq.gz,*_fastp.log}'
  publishDir "${params.outdir}/${sampleid}/03_fastqc_trimmed", mode: 'copy', pattern: '{*html,*json}'
  label "setting_4"

  input:
    tuple val(sampleid), path(fastq1), path(fastq2)

  output:
    path("${sampleid}.fastp.html")
    path("${sampleid}.fastp.json")
    path("${sampleid}_fastp.log")
    tuple val(sampleid), path("${sampleid}_1_qtrimmed.fastq.gz"), path("${sampleid}_2_qtrimmed.fastq.gz"), emit: trimmed_fq
    path("${sampleid}.fastp.json"), emit: fastp_json

  script:
    """
    fastp \
    --in1 ${fastq1} \
    --in2 ${fastq2} \
    --out1 ${sampleid}_1_qtrimmed.fastq.gz \
    --out2 ${sampleid}_2_qtrimmed.fastq.gz \
    --cut_front \
    --cut_tail \
    --json ${sampleid}.fastp.json \
    --html ${sampleid}.fastp.html \
    --thread 6 \
    --detect_adapter_for_pe \
    --length_required 50 --average_qual 20
    2>&1 | tee ${sampleid}_fastp.log
    """
}
//update to 2>&1 | tee ${sampleid}_fastp.log

process COVSTATS {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/05_mapping_to_consensus", mode: 'copy'

  input:
    tuple val(sampleid), path(bed), path(consensus), path(coverage), path(mapping_qual), path(top_hits), path(nanostats), val(target_size), path(reads_fasta), path(contig_seqids)
  output:
    path("*top_blast_with_cov_stats.txt")
    tuple val(sampleid), path("*top_blast_with_cov_stats.txt"), emit: detections_summary
    path("*top_blast_with_cov_stats.txt"), emit: detections_summary2

  script:
    """
    derive_coverage_stats.py --sample ${sampleid} --blastn_results ${top_hits} --nanostat ${nanostats} --coverage ${coverage} --bed ${bed} --target_size ${target_size} --contig_seqids ${contig_seqids} --reads_fasta ${reads_fasta} --consensus ${consensus} --mapping_quality ${mapping_qual}
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
    file("${sampleid}_megablast_top_viral_hits_filtered_with_contigs.txt")
    tuple val(sampleid), file("${sampleid}_megablast_top_viral_hits_filtered_with_contigs.txt"), emit: blast_results

  script:
    """
    fasta2table.py --fasta ${fasta} --sample ${sampleid} --tophits ${tophits}
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

process QCREPORT {
  publishDir "${params.outdir}/08_report", mode: 'copy', overwrite: true
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
  publishDir "${params.outdir}/logs", mode: 'copy', overwrite: true
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
  publishDir "${params.outdir}/${sampleid}/08_report", mode: 'copy', overwrite: true
  containerOptions "${bindOptions}"
  label 'setting_3'

  input:
    tuple val(sampleid), path(raw_nanoplot), path(filtered_nanoplot), path (rattle_status), path(consensus_fasta), path(top_blast_hits), path(blast_status), path(consensus_match_fasta), path(aln_sorted_bam), path(aln_sorted_bam_bai), path(blast_with_cov_stats),
    path(timestamp),
    path(qcreport_html),
    path(qcreport_txt),
    path(configyaml),
    path(samplesheet)

  output:
    path("*"), optional: true
    path("run_qc_report.html"), optional: true

  script:
  analyst_name = params.analyst_name.replaceAll(/ /, '_')
  facility = params.facility.replaceAll(/ /, '_')
    """
    cp ${qcreport_html} run_qc_report.html
    cp ${params.tool_versions} versions.yml
    cp ${params.default_params} default_params.yml

    build_report.py --samplesheet ${samplesheet} --result_dir . --params_file ${configyaml} --analyst ${analyst_name} --facility ${facility} --versions versions.yml --default_params_file default_params.yml
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

process EXTRACT_VIRAL_BLAST_HITS {
  tag "${sampleid}"
  label "setting_2"
  publishDir "${params.outdir}/${sampleid}/07_annotation", mode: 'copy'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(blast_results) 

  output:
    tuple val(sampleid), path("${sampleid}_megablast_top_viral_hits_filtered.txt"), emit: viral_blast_results
    path("${sampleid}_megablast_top_viral_hits.txt")
    path("${sampleid}_megablast_top_viral_hits_filtered.txt")
    path("${sampleid}_blastn.txt")

  script:
  """
  cat ${blast_results} > ${sampleid}_blastn.txt
  filter_blast.py --blastn_results ${sampleid}_blastn.txt --sample_name ${sampleid} --taxonkit_database_dir ${params.taxdump} --filter ${params.filter_terms}
  """
}

process EXTRACT_REF_FASTA {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/08_mapping_to_ref", mode: 'copy', pattern: '*fasta'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(blast_results)

  output:
    path("*fasta"), optional: true
    tuple val(sampleid), path("${sampleid}_ref_sequences.fasta"), emit: fasta_files, optional: true
  
  script:
    """
    cut -f1,3 ${blast_results} | sed '1d' | sed 's/ /_/g' > ids_to_retrieve.txt
    if [ -s ids_to_retrieve.txt ]
      then
        for i in `cut -f2  ids_to_retrieve.txt`; do j=`grep \${i} ids_to_retrieve.txt | cut -f1`; efetch -db nucleotide  -id \${i} -format fasta > ${sampleid}_\${i}.fasta ; done
        cat ${sampleid}_*.fasta > ${sampleid}_ref_sequences.fasta
    fi
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

process SAMTOOLS2 {
  publishDir "${params.outdir}/${sampleid}/08_mapping_to_ref", mode: 'copy'
  tag "${sampleid}"
  label 'setting_21'

  input:
    tuple val(sampleid), path(ref), path(sam)

  output:
    path "${sampleid}_ref_aln.sorted.bam"
    path "${sampleid}_ref_aln.sorted.bam.bai"
    path "${sampleid}_samtools_consensus_from_ref.fasta"
    path "${sampleid}_coverage.txt"
    tuple val(sampleid), path(ref), path("${sampleid}_ref_aln.sorted.bam"), path("${sampleid}_ref_aln.sorted.bam.bai"), emit: sorted_bam
    tuple val(sampleid), path("${sampleid}_coverage.txt"), emit: coverage
    tuple val(sampleid), path("${sampleid}_mapq.txt"), emit: mapping_quality

  script:
    """
    #filter out unmapped reads -F 4
    samtools view -@ ${task.cpus} -Sb -F 4 ${sam} | samtools sort -@ ${task.cpus} -o ${sampleid}_ref_aln.sorted.bam
    samtools index ${sampleid}_ref_aln.sorted.bam
    samtools coverage ${sampleid}_ref_aln.sorted.bam  > ${sampleid}_coverage.txt
    samtools view -@ ${task.cpus} ${sampleid}_ref_aln.sorted.bam | awk '{mapq[\$3]+=\$5; count[\$3]++} END {for (chr in mapq) printf "%s\\t%.2f\\n", chr, mapq[chr]/count[chr]}' > ${sampleid}_mapq.txt
    samtools consensus -@ ${task.cpus} -f fastq -a -A ${sampleid}_ref_aln.sorted.bam -o ${sampleid}_consensus.fastq 
    samtools consensus -@ ${task.cpus} -f fasta -a -A ${sampleid}_ref_aln.sorted.bam -o ${sampleid}_samtools_consensus_from_ref.fasta
    """
}

process BCFTOOLS {
  publishDir "${params.outdir}/${sampleid}/08_mapping_to_ref", mode: 'copy'
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
  publishDir "${params.outdir}/${sampleid}/08_mapping_to_ref", mode: 'copy'

  input:
   tuple val(sampleid), path(ref), path(bam), path(bai), path(vcf_applied_fasta)

  output:
    path("${sampleid}_bcftools_masked_consensus.fasta")
    path("${sampleid}_bcftools_masked_consensus.fasta"), emit: bcftools_masked_consensus_fasta

  script:
    """
		bedtools genomecov -ibam ${bam} -bga > ${sampleid}_genomecov.bed
    awk '\$4==0 {print}' ${sampleid}_genomecov.bed > ${sampleid}_zero_coverage.bed
    bedtools maskfasta -fi ${vcf_applied_fasta} -bed ${sampleid}_zero_coverage.bed -fo ${sampleid}_bcftools_masked_consensus.fasta
    """
}

process FASTQC_RAW {
  tag "$sampleid"
  publishDir "${params.outdir}/${sampleid}/01_fastqc_raw", mode: 'copy'
  label "setting_10"

  input:
    tuple val(sampleid), path(fastq1), path(fastq2)

  output:
    path("*_fastqc.{zip,html}")

  script:
    """
    fastqc --quiet --threads ${task.cpus} ${fastq1} ${fastq2}
    """
}

process FASTQC_TRIMMED {
    tag "$sampleid"
    label "setting_10"
    publishDir "${params.outdir}/${sampleid}/03_fastqc_trimmed", mode: 'copy'

    input:
      tuple val(sampleid), path(fastq1), path(fastq2)
    
    output:
      path("*_fastqc.{zip,html}")

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${fastq1} ${fastq2}
    """
}

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
    tuple val(sampleid), path("${sampleid}_non_rRNA_fwd.fastq.gz"), path("${sampleid}_non_rRNA_rev.fastq.gz"), emit: bbduk_filtered_fq

  script:
  """
  bbduk.sh -Xmx10g in=${fastq1} \
                   in2=${fastq2} \
                   out=${sampleid}_non_rRNA_fwd.fastq.gz \
                   out2=${sampleid}_non_rRNA_rev.fastq.gz \
                   outm=${sampleid}_rRNA_fwd.fastq.gz \
                   outm2=${sampleid}_rRNA_rev.fastq.gz \
                   k=31 ref=${ref} 2>${sampleid}_rRNA_reads.log
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
    path("${sampleid}_bbsplit.log")
    path("${sampleid}_cleaned_fwd.fastq.gz")
    path("${sampleid}_cleaned_rev.fastq.gz")
    path("${sampleid}_phyX_removed.fastq.gz")
    tuple val(sampleid), path("${sampleid}_cleaned_fwd.fastq.gz"), path("${sampleid}_cleaned_rev.fastq.gz"), emit: bbsplit_filtered_fq

  script:
  """
  bbsplit.sh -Xmx10g ref=${params.phix} \
             in=${fastq1} \
             in2=${fastq2} \
             out_phi-X174=${sampleid}_phyX_removed.fastq.gz \
             outu=${sampleid}_cleaned_fwd.fastq.gz \
             outu2=${sampleid}_cleaned_rev.fastq.gz 2> ${sampleid}_bbsplit.log
  """
}

process KRAKEN2 {
  tag "${sampleid}"
  label "setting_23"
  containerOptions "${bindOptions}"
  publishDir "${params.outdir}/${sampleid}/05_read_classification", mode: 'copy'

  input:
    tuple val(sampleid), path(fastq1), path(fastq2)
  output:
    path("${sampleid}_kraken2.log")
    path("${sampleid}_kraken2_report.txt")
    path("${sampleid}_kraken2_output.txt")
    path("${sampleid}_unclassified_1.fastq")
    path("${sampleid}_unclassified_2.fastq")
    tuple val(sampleid), path("${sampleid}_kraken2_report.txt"), path("${sampleid}_kraken2_output.txt"), path(fastq1), path(fastq2), path("${sampleid}_unclassified_1.fastq"), path("${sampleid}_unclassified_2.fastq"), emit: kraken2_results
    tuple val(sampleid), path("${sampleid}_kraken2_report.txt"), emit: kraken2_results2
  script:
  """
  kraken2 --db ${params.kraken2_db} --use-names \
          --paired --threads ${task.cpus} \
          --gzip-compressed \
          --confidence 0.05 \
          --report ${sampleid}_kraken2_report.txt \
          --output ${sampleid}_kraken2_output.txt \
          --unclassified-out ${sampleid}_unclassified#.fastq \
          --report-minimizer-data \
          --minimum-hit-groups 3 \
          ${fastq1} ${fastq2} > ${sampleid}_kraken2.log
  
  """
}

process BRACKEN {
  tag "${sampleid}"
  label 'setting_2'
  publishDir "$params.outdir/$sampleid/05_read_classification",  mode: 'copy'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(kraken_report)

  output:
    file("${sampleid}_bracken_report*.txt")

    tuple val(sampleid), path("${sampleid}_bracken_report_viral.txt"), emit: bracken_results

  script:
  """
  c1grep() { grep "\$@" || test \$? = 1; }

  est_abundance.py -i ${kraken_report} \
                  -k ${params.kraken2_db}/database50mers.kmer_distrib \
                  -t 1 \
                  -l S1 -o ${sampleid}_bracken_report.txt


  c1grep  "taxonomy_id\\|virus\\|viroid" ${sampleid}_bracken_report.txt > ${sampleid}_bracken_report_viral.txt
  awk -F'\\t'  '\$7>=0.0001'  ${sampleid}_bracken_report_viral.txt > ${sampleid}_bracken_report_viral_filtered.txt
  """
}

process KAIJU {
  tag "${sampleid}"
  label "setting_23"
  containerOptions "${bindOptions}"
  publishDir "${params.outdir}/${sampleid}/05_read_classification", mode: 'copy'

  input:
    tuple val(sampleid), path(fastq1), path(fastq2)
  output:
    file "${sampleid}_kaiju_name.tsv"
    file "${sampleid}_kaiju_summary*.tsv"
    file "${sampleid}_kaiju.krona"
    tuple val(sampleid), path("${sampleid}_kaiju_summary_viral.tsv"), emit: kaiju_results
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
  kaiju2table -e -t ${params.kaiju_nodes} -r species -n ${params.kaiju_names} -o ${sampleid}_kaiju_summary.tsv ${sampleid}_kaiju.tsv
  kaiju2krona -t ${params.kaiju_nodes} -n ${params.kaiju_names} -i ${sampleid}_kaiju.tsv -o ${sampleid}_kaiju.krona

  c1grep "taxon_id\\|virus\\|viroid\\|viricota\\|viridae\\|viriform\\|virales\\|virinae\\|viricetes\\|virae\\|viral" ${sampleid}_kaiju_summary.tsv > ${sampleid}_kaiju_summary_viral.tsv
  awk -F'\\t' '\$2>=0.05' ${sampleid}_kaiju_summary_viral.tsv > ${sampleid}_kaiju_summary_viral_filtered.tsv
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
    tuple val(sampleid), path("${sampleid}_cand_path_R1.fastq"), path("${sampleid}_cand_path_R2.fastq"), emit: fastq

  script:
  """
  extract_kraken_reads.py -k ${kraken_output} -r ${kraken_report} \
                          -t 10239 --include-children \
                          -s1 ${fastq1} -s2 ${fastq2} \
                          --fastq-output \
                          -o ${sampleid}_extracted_reads1.fq -o2 ${sampleid}_extracted_reads2.fq
  cat ${unc_fastq1} ${sampleid}_extracted_reads1.fq > ${sampleid}_cand_path_R1.fastq
  cat ${unc_fastq2} ${sampleid}_extracted_reads2.fq >  ${sampleid}_cand_path_R2.fastq
  """
}

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

process GENOMAD {
    tag "${sampleid}"
    label "setting_27"
    containerOptions "${bindOptions}"
    publishDir "${params.outdir}/${sampleid}/07_annotation", mode: 'copy'

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
      tuple val(sampleid), path("*_summary/*_virus.fna"), emit: virus_fasta
      

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
    """
}

/*
process RETRIEVE_VIRAL_READS_KRAKEN2 {
  tag "${sampleid}"
  label "setting_10"
  containerOptions "${bindOptions}"
  publishDir "${params.outdir}/${sampleid}/05_read_classification", mode: 'copy'

  input:
    tuple val(sampleid), path(kraken_report), path(kraken_output), path(fastq1), path(fastq2), path(unc_fastq1), path(unc_fastq2)
  output:
    tuple val(sampleid), path("${sampleid}_cand_path_R1.fastq"), path("${sampleid}_cand_path_R2.fastq"), emit: fastq

  script:
  """
  extract_kraken_reads.py -k ${kraken_output} -r ${kraken_report} \
                          -t 10239 --include-children \
                          -s1 ${fastq1} -s2 ${fastq2} \
                          --fastq-output \
                          -o ${sampleid}_extracted_reads1.fq -o2 ${sampleid}_extracted_reads2.fq
  cat ${unc_fastq1} ${sampleid}_extracted_reads1.fq > ${sampleid}_cand_path_R1.fastq
  cat ${unc_fastq2} ${sampleid}_extracted_reads2.fq >  ${sampleid}_cand_path_R2.fastq
  """
}
/*
process MERGE_WITH_UNCLASSIFIED_READS_KRAKEN2 {
  tag "${sampleid}"
  label "setting_10"
  containerOptions "${bindOptions}"
  publishDir "${params.outdir}/${sampleid}/05_read_classification", mode: 'copy'

  input:
    tuple val(sampleid), path(kraken_report), path(kraken_output), path(fastq1), path(fastq2), path(unc_fastq1), path(unc_fastq2)
  output:
        tuple val(sampleid), path("${sampleid}_cand_path_R1.fastq"), path("${sampleid}_cand_path_R2.fastq"), emit: fastq

  script:
  """
  grep 

  """
}
*/


workflow {
  TIMESTAMP_START ()
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

  FASTP ( ch_sample )
  FASTQC_RAW (ch_sample )
  FASTQC_TRIMMED (FASTP.out.trimmed_fq )
  trial_ch = FASTP.out.trimmed_fq.map { sample_id, read1, read2 ->
    tuple(sample_id, read1, read2, file(params.sortmerna_ref), file(params.sortmerna_idx))
  }
  trial_ch2 = FASTP.out.trimmed_fq.map { sample_id, read1, read2 ->
    tuple(sample_id, read1, read2, file(params.sortmerna_ref))
  }

  //Filtering with sortmerna takes much longer than bbduk
  //SORTMERNA ( trial_ch )
  BBDUK ( trial_ch2 )
  //remove phiX reads
  FILTER_CONTROL ( BBDUK.out.bbduk_filtered_fq )
  //provide option to filter host
  //filter host by default?
  // read classification with Kraken
  KRAKEN2 ( FILTER_CONTROL.out.bbsplit_filtered_fq )
  BRACKEN ( KRAKEN2.out.kraken2_results2 )
  //retrieve reads that were not classified and reads classified as viral
  RETRIEVE_VIRAL_READS_KRAKEN2 ( KRAKEN2.out.kraken2_results )
  //merge unclassified reads with viral reads from kraken2
  //MERGE_WITH_UNCLASSIFIED_READS_KRAKEN2 ( RETRIEVE_VIRAL_READS_KRAKEN2.out.fastq )

  KAIJU ( FILTER_CONTROL.out.bbsplit_filtered_fq )
  
  //perform de novo assembly with spades using rnaspades
  SPADES ( RETRIEVE_VIRAL_READS_KRAKEN2.out.fastq )
   


  //SPADES ( FILTER_CONTROL.out.bbsplit_filtered_fq )
  //SPADES ( SORTMERNA.out.fastp_filtered_fq )
  SEQTK ( SPADES.out.assembly )
  //DIAMOND  ( SEQTK.out.filt_fasta.splitFasta(by: 5000, file: true) )
  //DIAMOND.out.diamond_results
  //  .groupTuple()
  //  .set { ch_blastxresults } 
  //At the moment we are not 
  BLASTN( SEQTK.out.filt_fasta.splitFasta(by: 5000, file: true) )
  BLASTN.out.blast_results
    .groupTuple()
    .set { ch_blastresults } 
  EXTRACT_VIRAL_BLAST_HITS ( ch_blastresults )
  //Add consensus sequence to blast results summary table
  FASTA2TABLE ( EXTRACT_VIRAL_BLAST_HITS.out.viral_blast_results.join(SEQTK.out.filt_fasta) )
  EXTRACT_REF_FASTA ( EXTRACT_VIRAL_BLAST_HITS.out.viral_blast_results )
  mapping_ch = EXTRACT_REF_FASTA.out.fasta_files.join(FILTER_CONTROL.out.bbsplit_filtered_fq)
  MAPPING_BACK_TO_REF ( mapping_ch )
  SAMTOOLS2 ( MAPPING_BACK_TO_REF.out.aligned_sam )
  BCFTOOLS ( SAMTOOLS2.out.sorted_bam )
  BEDTOOLS ( BCFTOOLS.out.vcf_applied_fasta )
  PYFAIDX ( EXTRACT_REF_FASTA.out.fasta_files )
  MOSDEPTH (SAMTOOLS2.out.sorted_bam.join(PYFAIDX.out.bed))
  
  cov_stats_summary_ch = MOSDEPTH.out.mosdepth_results.join(FASTA2TABLE.out.blast_results)
                                                      .join(FASTP.out.fastp_json)
                                                      .join(BEDTOOLS.out.bcftools_masked_consensus_fasta)
                                                      .join(SAMTOOLS2.out.coverage)
                                                      .join(SAMTOOLS2.out.mapping_quality)
                                                      .join(EXTRACT_REF_FASTA.out.fasta_files)
                                                      
  COVSTATS( cov_stats_summary_ch )
  
  GENOMAD ( SPADES.out.assembly )


  
  //MOSDEPTH (SAMTOOLS2.out.sorted_bams.join(PYFAIDX.out.bed))
  //trimmed_fq = FASTP.out.fastp_trimmed_fq

    // Perform quality filtering of reads using chopper
//    if (params.qual_filt) {
//      CHOPPER ( trimmed_fq)f
//      filtered_fq = CHOPPER.out.chopper_filtered_fq
//    }
//    else { filtered_fq = trimmed_fq
//    }

    //Reformat fastq read names after the first whitespace
//    REFORMAT( filtered_fq )
/*
    //Run Nanoplot on merged raw fastq files after data processing
  QC_POST_DATA_PROCESSING ( FASTP.out.fastp_trimmed_fq )

    if (params.subsample) {
      SUBSAMPLE ( REFORMAT.out.reformatted_fq )
      final_fq = SUBSAMPLE.out.subsampled_fq
    }
    else {
      final_fq = REFORMAT.out.reformatted_fq
    }

    //Derive QC report
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(QC_PRE_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_POST_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
    QCREPORT(ch_multiqc_files.collect())

    if (!params.preprocessing_only) {
      //Currently only one analysis mode in ont_amplicon, consider removing if no other mode is added to this pipeline
      //We have had talks about including an option to just a map to a reference of interest, but this is not implemented yet
      if ( params.analysis_mode == 'clustering' ) {
        //Perform clustering using Rattle and convert to fasta file
        //Branch outputs of Rattle into passed and failed
        //If the clustering step succeeds, it will proceed to the polishing step
        ch_fq_target_size = (final_fq.join(ch_target_size))
        RATTLE ( ch_fq_target_size )
        CLUSTER2FASTA ( RATTLE.out.clusters )

        //Polish consensus sequence using Racon followed by Medaka and samtools consensus
        if (params.polishing) {
          MINIMAP2_RACON ( CLUSTER2FASTA.out.fasta )
          RACON ( MINIMAP2_RACON.out.draft_mapping)
          MEDAKA2 ( RACON.out.polished )
          ch_branched = MEDAKA2.out.consensus2
            | branch { sampleid, clusters, consensus, status ->
              passed: status == "passed"
              failed: status == "failed"
            }

          ch_passed = ch_branched.passed
            | map { sampleid, clusters, consensus, status -> [sampleid, consensus] }
          ch_failed = ch_branched.failed
            | map { sampleid, clusters, consensus, status -> [sampleid, clusters] }

          consensus = ch_passed.concat(ch_failed.ifEmpty([]))
        }
        //If polishing is skipped, directly use the clusters generated by Rattle for blast search
        else {
          consensus = CLUSTER2FASTA.out.fasta2
        }

        //Remove trailing Ns and primer sequences from consensus sequence
          CUTADAPT ( consensus.join(ch_primers) )

        //Blast steps for samples targetting COI
        ch_coi_for_blast = (CUTADAPT.out.trimmed.join(ch_coi))
        //Blast to COI database
        BLASTN_COI(ch_coi_for_blast)
        //Identify consensus that are in the wrong orientation and reverse complement them
        ch_revcomp = (CUTADAPT.out.trimmed.join(BLASTN_COI.out.coi_blast_results))
        REVCOMP ( ch_revcomp )
        //Blast to NCBI nt database
        BLASTN ( REVCOMP.out.revcomp )

        //Directly blast to NCBI nt database all other samples
        ch_other_for_blast = (CUTADAPT.out.trimmed.join(ch_other))
        BLASTN2 ( ch_other_for_blast )

        //Merge blast results from all samples
        ch_blast_merged = BLASTN.out.blast_results.mix(BLASTN2.out.blast_results.ifEmpty([]))

        ch_blast_merged2 = ch_blast_merged.map { sampleid, blast_results, status -> [sampleid, blast_results] }

        //Extract top blast hit, assign taxonomy information to identify consensus that match target organism
        EXTRACT_BLAST_HITS ( ch_blast_merged2.join(ch_targets) )
        //Add consensus sequence to blast results summary table
        FASTA2TABLE ( EXTRACT_BLAST_HITS.out.topblast.join(consensus) )

        //MAPPING BACK TO CONSENSUS
        mapping2consensus_ch = (EXTRACT_BLAST_HITS.out.consensus_fasta_files.join(REFORMAT.out.cov_derivation_ch))
        //Map filtered reads back to the portion of sequence which returned a blast hit
        MINIMAP2_CONSENSUS ( mapping2consensus_ch )
        //Derive bam file and coverage statistics
        SAMTOOLS_CONSENSUS ( MINIMAP2_CONSENSUS.out.aligned_sample )
        //Derive bed file for mosdepth to run coverage statistics
        PYFAIDX ( EXTRACT_BLAST_HITS.out.consensus_fasta_files )
        MOSDEPTH (SAMTOOLS_CONSENSUS.out.sorted_bams.join(PYFAIDX.out.bed))
        SEQTK (SAMTOOLS_CONSENSUS.out.contig_seqids.join(final_fq))
        //Derive summary file presenting coverage statistics alongside blast results
        cov_stats_summary_ch = MOSDEPTH.out.mosdepth_results.join(EXTRACT_BLAST_HITS.out.consensus_fasta_files)
                                                             .join(SAMTOOLS_CONSENSUS.out.coverage)
                                                             .join(SAMTOOLS_CONSENSUS.out.mapping_quality)
                                                             .join(FASTA2TABLE.out.blast_results)
                                                             .join(QC_POST_DATA_PROCESSING.out.filtstats)
                                                             .join(ch_target_size)
                                                             .join(SEQTK.out.fasta)
                                                             .join(SEQTK.out.contig_seqids)

        COVSTATS(cov_stats_summary_ch)

        files_for_report_ind_samples_ch = QC_PRE_DATA_PROCESSING.out.rawnanoplot.join((QC_POST_DATA_PROCESSING.out.filtnanoplot)
                                                                                .join(RATTLE.out.status)
                                                                                .join(CUTADAPT.out.trimmed)
                                                                                .join(ch_blast_merged)
                                                                                .join(SAMTOOLS_CONSENSUS.out.sorted_bams)
                                                                                .join(COVSTATS.out.detections_summary))
        files_for_report_global_ch = TIMESTAMP_START.out.timestamp
            .concat(QCREPORT.out.qc_report_html)
            .concat(QCREPORT.out.qc_report_txt)
            .concat(configyaml)
            .concat(Channel.from(params.samplesheet).map { file(it) }).toList()

        HTML_REPORT(files_for_report_ind_samples_ch
            .combine(files_for_report_global_ch))

        //MAPPING BACK TO REFERENCE
        if (params.mapping_back_to_ref) {
          mapping_ch = (EXTRACT_BLAST_HITS.out.reference_fasta_files.join(REFORMAT.out.cov_derivation_ch))
          //Map filtered reads back to the reference sequence which was retrieved from blast search
          MINIMAP2_REF ( mapping_ch )
          //Derive bam file and consensus fasta file
          SAMTOOLS ( MINIMAP2_REF.out.aligned_sample )
        }
      }
/*
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
