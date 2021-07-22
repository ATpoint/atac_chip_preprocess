// Trim fastq files with cutadapt and stream into bowtie2 for alignment:

process Bowtie2Align {

    tag "$sample_id"

    cpus params.threads
    memory params.memory

    publishDir params.align_dir, mode: params.align_pubmode

    input:
    tuple val(sample_id), path(reads)
    path(idx)                         

    output:
    tuple val("${sample_id}"), path("${sample_id}_raw.bam"), path("${sample_id}_raw.bam.bai"), emit: bam

    path("${sample_id}_raw.flagstat"), emit: flagstat
    
    script:
    
    """

    if [[ $params.idx == '' ]]; then
        use_idx=`echo $idx | awk -F \".\" '{print \$1 | \"sort -u\"}'`
    else use_idx=$params.idx
    fi    

    bowtie2_constant='bowtie2 -q --threads $params.align_threads --rg-id $sample_id $params.align_additional -x \$use_idx'

    trim_constant='cutadapt $params.trim_additional --quiet -j 1 -a $params.trim_adapter'

    sort_constant='samtools sort $params.sort_additional --write-index -@ $params.sort_threads -m $params.sort_mem -o ${sample_id}_raw.bam##idx##${sample_id}_raw.bam.bai'

    if [[ $params.mode == "paired" ]]; then

        eval \$trim_constant -A $params.trim_adapter --interleaved -o - ${reads[0]} ${reads[1]} \
        | eval \$bowtie2_constant --interleaved - \
        | samblaster --ignoreUnmated \
        | eval \$sort_constant

    fi

    if [[ $params.mode == "single" ]]; then

        \$trim_constant ${reads} \
        | eval \$bowtie2_constant -U - \
        | samblaster --ignoreUnmated \
        | eval \$sort_constant

    fi

    samtools flagstat -@ $params.align_threads ${sample_id}_raw.bam > ${sample_id}_raw.flagstat
    
    """

} 