#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { VIEW                  } from './workflows/view/main'


def helpMessage () {
    log.info """
    Virus Integrated Evaluation Workflow
    Marie-Emilie Gauthier
    Cameron Hyde

    Usage:
    Run the command
    nextflow run main.nf -profile singularity -params-file {params.yml}

    Required arguments:
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

workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }

    VIEW ()
}
