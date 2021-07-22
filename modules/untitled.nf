// Trim fastq files with cutadapt and stream into bowtie2 for alignment:

process Bowtie2Align {

    tag "$sample_id"

    cpus   params.threads
    memory params.mem

    publishDir params.outdir, mode: params.pubmode

    input:
    tuple val(sample_id), path(reads)
    path(idx)                         

    output:
    path("${sample_id}_raw.bam")
    
    script:
    // inspired from https://github.com/nf-core/modules/blob/master/modules/bowtie2/align/main.nf
    """

    INDEX=`find -L $idx -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`

    cut_basic=\"cutadapt --quiet -j 1 -m 18 --max-n 0.1 -a\"
    #/ case: paired-end data with trmimming:
    if [[ $params.mode == "paired" ]]; then
        seqtk mergepe ${reads[0]} -2 ${reads[1]} \
        | eval \$cut_basic -a $params.adapter -A $params.adapter $params.trim_additional --interleaved - \
        | bowtie2 -q --threads $params.threads --rg-id $sample_id $params.align_additional -x \$INDEX --interleaved - \
        | samtools view -o \"${sample_id}_raw.bam\"
    fi

    if [[ $params.mode == "single" ]]; then
        eval \$cut_basic -a $params.adapter $params.trim_additional $reads \
        | bowtie2 -q --threads $params.threads --rg-id $sample_id $params.align_additional -x \$INDEX - \
        | samtools view -o \"${sample_id}_raw.bam\"
    fi    

    """

} 

#| samblaster --quiet --ignoreUnmated - - \
    #| samtools sort $params.sort_additional --write-index -m $params.sort_mem -@ $params.sort_threads -o ${sample_id}_raw.bam         