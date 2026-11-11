//might want to specify parameter k=31 outside of process in the future
process BBMAP {
    tag "$meta.id"
    label 'setting_11'
    publishDir { "${params.outdir}/${meta.id}/04_cleaned" }, mode: 'copy', pattern: '{*cleaned*fastq.gz,*bbsplit.log,*bbsplit_stats.txt,*bbduk.log}'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"
    
    input:
    tuple val(meta), path(reads)
    val(mounted_db)
    path(staged_db)
    val(subsample_enabled)
    val(sample_size)
    path  index, name: 'input_index'
    path  primary_ref
    tuple val(other_ref_names), path(other_ref_paths)
    val   only_build_index

    output:
    //path("${meta.id}_non_rRNA_1.fastq.gz")
    //path("${meta.id}_non_rRNA_2.fastq.gz")
    path("${task.ext.prefix ?: meta.id}_bbduk.log")
    path("${task.ext.prefix ?: meta.id}_bbsplit.log"), optional: true
    path("${task.ext.prefix ?: meta.id}_subsampling.log"), optional: true
    
    //path('*_bbduk.log')                      , emit: log2
    path('*_bbduk.log')            , emit: bbduk_log
    path "bbsplit_index"                      , optional:true, emit: index
    tuple val(meta), path('*primary*fastq.gz'), optional:true, emit: primary_fastq
    //tuple val(meta), path('*fastq.gz')        , emit: all_fastq
    tuple val(meta), path('*subsampled*.fastq.gz') , emit: all_fastq
    path('*bbsplit_stats.txt')             , emit: stats
    //path('*bbsplit_stats.txt')                , emit: stats2
    //tuple val(meta), path('*.log')            , optional:true, emit: log
    tuple val(meta), path('*_analysed_read_count.txt'), emit: analysed_read_count
    
    path "versions.yml"                       , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def rRNA_clean_out  = meta.single_end ? "out=${prefix}_non_rRNA.fastq.gz" : "out1=${prefix}_non_rRNA_1.fastq.gz out2=${prefix}_non_rRNA_2.fastq.gz"
    def rRNA_clean_in = meta.single_end ? "in=${prefix}_non_rRNA.fastq.gz" : "in1=${prefix}_non_rRNA_1.fastq.gz in2=${prefix}_non_rRNA_2.fastq.gz"
    def subsampled = meta.single_end ? "out=${prefix}_subsampled.fastq.gz" : "out1=${prefix}_subsampled_1.fastq.gz out2=${prefix}_subsampled_2.fastq.gz"
    def db = mounted_db ?: staged_db
    def contaminants_fa = db ? "ref=${db}" : ''
    if ( subsample_enabled && !sample_size ) {
        error "The bbmap_bbduk process must have a sample_size value included"
    }

    def other_refs = []
    other_ref_names.eachWithIndex { name, i ->
        other_refs << "ref_${name}=${other_ref_paths[i]}"
    }
    def avail_mem = (task.memory.toMega() * 0.85).toInteger()
    def index_files=''
    def refstats_cmd=''
    def use_index = index ? true : false

    def phiX_clean_out = meta.single_end
            ? "outu=${prefix}_cleaned.fastq.gz"
            : "outu=${prefix}_cleaned_1.fastq.gz outu2=${prefix}_cleaned_2.fastq.gz"

    def phiX_clean_in = meta.single_end
        ? "in=${prefix}_cleaned.fastq.gz"
        : "in1=${prefix}_cleaned_1.fastq.gz in2=${prefix}_cleaned_2.fastq.gz"

    def primary_out = "out_primary=${prefix}_phiX_aligning.fastq.gz"
    refstats_cmd = 'statsfile=' + prefix + '_bbsplit_stats.txt'

    if (only_build_index) {
        if (primary_ref && other_ref_names && other_ref_paths) {
            index_files = 'ref_primary=' +primary_ref + ' ' + other_refs.join(' ') + ' path=bbsplit_build'
        } else {
            log.error 'ERROR: Please specify as input a primary fasta file along with names and paths to non-primary fasta files.'
        }
    } else {
        if (index) {
            index_files = "path=index_writable"
        //} else if (primary_ref && other_ref_names && other_ref_paths) {
        } else if (primary_ref) {
            def other_refs_str = other_ref_names ? other_refs.join(' ') : ''
            index_files = "ref_primary=${primary_ref} ${other_refs_str}".trim()
        } else {
            log.error 'ERROR: Please specify either a BBSplit index or a primary fasta file.'
        }
        //fastq_out = meta.single_end ? "basename=${prefix}_%.fastq.gz" : "basename=${prefix}_%_#.fastq.gz"

        

        //refstats_cmd = 'refstats=' + prefix + '_bbsplit_stats.txt'
        //statsfile=${sampleid}_bbsplit_stats.txt
        
    }
    """
    # If using a pre-built index, create writable structure: symlink all files except
    # summary.txt (which we copy to modify). When we stage in the index files the time
    # stamps get disturbed, which bbsplit doesn't like. Fix the time stamps in summaries.
    if [ "$use_index" == "true" ]; then
        find input_index/ref -type f | while read -r f; do
            target="index_writable/\${f#input_index/}"
            mkdir -p "\$(dirname "\$target")"
            [[ \$(basename "\$f") == "summary.txt" ]] && cp "\$f" "\$target" || ln -s "\$(realpath "\$f")" "\$target"
        done
        find index_writable/ref/genome -name summary.txt | while read -r summary_file; do
            src=\$(grep '^source' "\$summary_file" | cut -f2- -d\$'\\t' | sed 's|.*/ref/|index_writable/ref/|')
            mod=\$(echo "System.out.println(java.nio.file.Files.getLastModifiedTime(java.nio.file.Paths.get(\\"\$src\\")).toMillis());" | jshell -J-Djdk.lang.Process.launchMechanism=vfork - 2>/dev/null | grep -oE '^[0-9]{12,14}\$')
            sed -e 's|bbsplit_index/ref|index_writable/ref|' -e "s|^last modified.*|last modified\\t\$mod|" "\$summary_file" > \${summary_file}.tmp && mv \${summary_file}.tmp \${summary_file}
        done
    fi

    # Filter rRNA
    bbduk.sh \\
        -Xmx${avail_mem}M \\
        $raw \\
        $rRNA_clean_out \\
        k=31 \\
        threads=${task.cpus} \\
        $args \\
        $contaminants_fa \\
        &>${prefix}_bbduk.log


    # Filter PhiX
    bbsplit.sh \\
        -Xmx${avail_mem}M \\
        $index_files \\
        threads=$task.cpus \\
        $rRNA_clean_in \\
        $phiX_clean_out \\
        $primary_out \\
        $refstats_cmd \\
        $args 2>| >(tee ${prefix}_bbsplit.log >&2)

    # Summary files will have an absolute path that will make the index
    # impossible to use in other processes - fix paths and rename automically
    if [ -d bbsplit_build/ref/genome ]; then
        find bbsplit_build/ref/genome -name summary.txt | while read -r summary_file; do
            sed "s|^source.*|source\\t\$(grep '^source' "\$summary_file" | cut -f2- -d\$'\\t' | sed 's|.*/bbsplit_build|bbsplit_index|')|" "\$summary_file" > \${summary_file}.tmp && mv \${summary_file}.tmp \${summary_file}
        done
        mv bbsplit_build bbsplit_index
    fi

    
    if [ "${subsample_enabled}" = "true" ]; then
        # Read count
        THRESHOLD=$sample_size
        if [ "${meta.single_end}" = "true" ]; then
            READS=\$(gzip -cd ${prefix}_cleaned.fastq.gz | awk 'END {print NR/4}')
        else
            READS=\$(gzip -cd ${prefix}_cleaned_1.fastq.gz | awk 'END {print NR/4}')
        fi
        

        #If read counts exceed the threshold, perform subsampling
        if [ "\$READS" -gt "\$THRESHOLD" ]; then
            SAMPLERATE=\$(awk -v target="\$THRESHOLD" -v total="\$READS" \
                'BEGIN { printf "%.10f", target / total }')
            
            bbduk.sh \\
            -Xmx${avail_mem}M \\
            $phiX_clean_in \\
            $subsampled \\
            samplerate=\$SAMPLERATE \\
            sampleseed=100 \\
            &>${prefix}_subsampling.log

            
        else
            if [ "${meta.single_end}" = "true" ]; then
                ln -s "${prefix}_cleaned.fastq.gz" "${prefix}_subsampled.fastq.gz"
            else
                ln -s "${prefix}_cleaned_1.fastq.gz" "${prefix}_subsampled_1.fastq.gz"
                ln -s "${prefix}_cleaned_2.fastq.gz" "${prefix}_subsampled_2.fastq.gz"
            fi
        fi
    else
        if [ "$meta.single_end" = true ]; then
            ln -s "${prefix}_cleaned.fastq.gz" "${prefix}_subsampled.fastq.gz"
        else
            ln -s "${prefix}_cleaned_1.fastq.gz" "${prefix}_subsampled_1.fastq.gz"
            ln -s "${prefix}_cleaned_2.fastq.gz" "${prefix}_subsampled_2.fastq.gz"
        fi
    fi

    if [ "${meta.single_end}" = "true" ]; then
        ANALYSED_READS=\$(gzip -cd ${prefix}_subsampled.fastq.gz | awk 'END {print NR/4}')
    else
        ANALYSED_READS=\$(gzip -cd ${prefix}_subsampled_1.fastq.gz | awk 'END {print NR/4}')
    fi

    echo "\$ANALYSED_READS" > ${prefix}_analysed_read_count.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "]")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_command  = meta.single_end ? "echo '' | gzip > ${prefix}.fastq.gz" : "echo '' | gzip > ${prefix}_1.fastq.gz ; echo '' | gzip > ${prefix}_2.fastq.gz"
    def other_refs = ''
    other_ref_names.eachWithIndex { name, i ->
        other_refs += "echo '' | gzip > ${prefix}_${name}.fastq.gz"
    }
    def will_build_index = only_build_index || (!index && primary_ref && other_ref_names && other_ref_paths)
    """
    touch ${prefix}_bbduk.log
    $output_command

    # Create index directory if building an index (either only_build_index or on-the-fly)
    if [ "${will_build_index}" == "true" ]; then
        mkdir -p bbsplit_index
    fi

    # Only create output files if splitting (not just building index)
    if ! (${only_build_index}); then
        echo '' | gzip >  ${prefix}_primary.fastq.gz
        ${other_refs}
        touch ${prefix}.stats.txt
    fi

    touch ${prefix}_bbsplit.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS

    """
}